clc;clear;

% ��֤���ġ���̬��Ҷ˹����������������ѡȡ�������ķ���ʵ�鲿�֡���linprog()����
% ������ʣ���������ð���͵ĸ���
p = 0.6;
P_A = [p,1-p];      % �����߶ԡ������ߡ����͵���������
P_D = [0.3,0.4,0.4];        % �����߶ԡ������ߡ����͵���������
num_type_A = 2;     % ������������
num_type_D = 1;     % ������������
num_type = [num_type_A, num_type_D];    % ��1�����ǹ���������������2�����Ƿ�����������

% ����ð���͹���������µ��������
high_benefit_matrix = {
    [-730,1050], [-680,980], [-580,950];
    [-750,1200], [-740,1000], [-680,1100];
    [-560,900], [-700,980], [-680,1000]
};

% ���屣���͹���������µ��������
low_benefit_matrix = {
    [-780,1280],[-500,810],[-540,900];
    [-700,1100],[-670,1100],[-610,1080];
    [-510,1000],[-580,1020],[-580,1180];
};

% ��������������ڼ���������������
U = cat(3, high_benefit_matrix, low_benefit_matrix);

% ��ȡ D �Ĳ��������������������
nD = size(high_benefit_matrix, 1);
% ��ȡ A �Ĳ��������������������
nA = size(high_benefit_matrix, 2);

% ��ʼ�����Ը���
F_A_case1 = [0,0.3417,0.6583];
F_A_case2 = [0.8547,0,0.1453];
F_star_A = [F_A_case1;F_A_case2];
F_star_D = [0,0,1];

% ��ʼ������������󣬷ֱ�洢 D �� A ����������
expected_D_matrix = zeros(nD, nA);
expected_A_matrix = zeros(nD, nA);

% ���������������
for t_A = 1:num_type_A  % �������п��ܵĹ��������� t_A
    for i = 1:nD
        for j = 1:nA
            expected_D_matrix(i, j) = p * high_benefit_matrix{i, j}(1) * F_star_D(i) + (1 - p) * low_benefit_matrix{i, j}(1) * F_star_D(i);
            expected_A_matrix(i, j) = p * high_benefit_matrix{i, j}(2) * F_star_A(t_A,j) + (1 - p) * low_benefit_matrix{i, j}(2) * F_star_A(t_A,j);
        end
    end
end

% ����������������
efficiency_D = zeros(1,nD);
efficiency_A = zeros(1,nA);

% ���¾����µĲ���
for k = 1:nD % �����߸������ԵĲ�������
    efficiency_D(k) = calculate_efficiency_ver2(num_type, P_A, F_star_A, U, k, 'defender');
end
for k = 1:nA % �����߸������ԵĲ�������
    efficiency_A(k) = calculate_efficiency_ver2(num_type, P_D, F_star_D, U, k, 'attacker');
end

% ���������������֤�������������ĵ���ȷ��
disp('�����߲��ԵĹ���Ч��:');
disp(efficiency_A);
disp('�����߲��Եķ���Ч��:');
disp(efficiency_D);

% ������Ʋ��Բ��� - D
dominant_D_strategy = false(nD, 1);
for i = 1:nD
    is_dominant = true;
    for j = 1:nD
        if i ~= j
            for k = 1:nA
                if expected_D_matrix(i, k) < expected_D_matrix(j, k)
                    is_dominant = false;
                    break;
                end
            end
            if ~is_dominant
                break;
            end
        end
    end
    dominant_D_strategy(i) = is_dominant;
end

% ������Ʋ��Բ��� - A
dominant_A_strategy = false(nA, 1);
for i = 1:nA
    is_dominant = true;
    for j = 1:nA
        if i ~= j
            for k = 1:nD
                if expected_A_matrix(k, i) < expected_A_matrix(k, j)
                    is_dominant = false;
                    break;
                end
            end
            if ~is_dominant
                break;
            end
        end
    end
    dominant_A_strategy(i) = is_dominant;
end

% ���û�����Ʋ��ԣ�ʹ�������������������
if ~any(dominant_D_strategy) && ~any(dominant_A_strategy)
    if strcmp(mode, 'mixed')
        % D ��ϲ������
        f_D = [-1, zeros(1, nD)];
        A_D = [ones(nD, 1), -expected_D_matrix];
        b_D = zeros(nD, 1);
        Aeq_D = [0, ones(1, nD)];
        beq_D = 1;
        lb_D = [-inf, zeros(1, nD)];
        ub_D = [inf, ones(1, nD)];
        options = optimoptions('linprog', 'Display', 'off');
        [x_D, ~, exitflag_D] = linprog(f_D, A_D, b_D, Aeq_D, beq_D, lb_D, ub_D, options);
        if exitflag_D ~= 1
            error('D ��ϲ������Թ滮���ʧ��');
        end
        F_star_D = x_D(2:end);
        D_max_expected_payoff = -x_D(1);

        % A ��ϲ������
        f_A = [-1, zeros(1, nA)];
        A_A = [ones(nD, 1), expected_A_matrix'];
        b_A = zeros(nD, 1);
        Aeq_A = [0, ones(1, nA)];
        beq_A = 1;
        lb_A = [-inf, zeros(1, nA)];
        ub_A = [inf, ones(1, nA)];
        [x_A, ~, exitflag_A] = linprog(f_A, A_A, b_A, Aeq_A, beq_A, lb_A, ub_A, options);
        if exitflag_A ~= 1
            error('A ��ϲ������Թ滮���ʧ��');
        end
        F_star_A = x_A(2:end);
        A_min_expected_payoff = -x_A(1);

        % ����������������
        efficiency_D = zeros(1,nD);
        efficiency_A = zeros(1,nA);
        for k = 1:nD
            efficiency_D(k) = calculate_efficiency_ver2(num_type, P_A, F_star_A, U, k, 'defender');
        end
        for k = 1:nA
            efficiency_A(k) = calculate_efficiency_ver2(num_type, P_D, F_star_D, U, k, 'attacker');
        end

        disp('D �Ļ�ϲ���:');
        disp(F_star_D);
        disp(['D �������������: ', num2str(D_max_expected_payoff)]);
        D_opt_strategy = find(F_star_D);
        disp('A �Ļ�ϲ���:');
        disp(F_star_A);
        disp(['A �������������: ', num2str(A_min_expected_payoff)]);
        A_opt_strategy = find(F_star_A);

        % ��������ȷ���������Դ�������
        A_strategies = cell(1, nA);
        D_strategies = cell(1, nD);
        for i = 1:nA
            A_strategies{i} = sprintf('a%d', i);
        end
        for i = 1:nD
            D_strategies{i} = sprintf('d%d', i);
        end
        optimal_A_strategy = A_strategies(A_opt_strategy);
        optimal_D_strategy = D_strategies(D_opt_strategy);

        disp(['A �������Ų���: ', strjoin(optimal_A_strategy, ', '),'����ϲ��Ը��ʣ�',num2str(F_star_A')]);
        disp(['A ���Ĺ�������Ч��������',num2str(efficiency_A)]);
        disp(['D �������Ų���: ', strjoin(optimal_D_strategy, ', '),'����ϲ��Ը��ʣ�',num2str(F_star_D')]);
        disp(['D ���Ĺ�������Ч��������',num2str(efficiency_D)]);
    elseif strcmp(mode, 'pure')
        % ��ʼ�� A �� D �����Ų�������
        A_opt_strategy = 1;
        D_opt_strategy = 1;

        % ��ʼ�������������
        max_A_income = -Inf;
        max_D_income = -Inf;

        % Ѱ�� A �����Ų��ԣ�����������棩
        for i = 1:nA
            total_income = sum(expected_A_matrix(:, i));
            if total_income > max_A_income
                max_A_income = total_income;
                A_opt_strategy = i;
            end
        end

        % Ѱ�� D �����Ų��ԣ�����������棩
        for j = 1:nD
            total_income = sum(expected_D_matrix(j, :));
            if total_income > max_D_income
                max_D_income = total_income;
                D_opt_strategy = j;
            end
        end

        % ����������������
        efficiency_D = zeros(1,nD);
        efficiency_A = zeros(1,nA);
        F_star_A = zeros(1, nA);
        F_star_A(A_opt_strategy) = 1;
        F_star_D = zeros(1, nD);
        F_star_D(D_opt_strategy) = 1;
        for k = 1:nD
            efficiency_D(k) = calculate_efficiency(num_type, P_A, F_star_A, U, k, 'defender');
        end
        for k = 1:nA
            efficiency_A(k) = calculate_efficiency(num_type, P_D, F_star_D, U, k, 'attacker');
        end

        % ��������ȷ���������Դ�������
        A_strategies = cell(1, nA);
        D_strategies = cell(1, nD);
        for i = 1:nA
            A_strategies{i} = sprintf('a%d', i);
        end
        for i = 1:nD
            D_strategies{i} = sprintf('d%d', i);
        end
        optimal_A_strategy = A_strategies{A_opt_strategy};
        optimal_D_strategy = D_strategies{D_opt_strategy};

        disp(['A �������Ų���: ', num2str(optimal_A_strategy)]);
        disp(['A ���Ĺ�������Ч��������',num2str(efficiency_A)]);
        disp(['D �������Ų���: ', num2str(optimal_D_strategy)]);
        disp(['D ���Ĺ�������Ч��������',num2str(efficiency_D)]);
    end
else
    % ��������Ʋ��ԣ���ֱ�Ӹ������Ʋ���ȷ�����Ų���
    if any(dominant_D_strategy)
        D_opt_strategy = find(dominant_D_strategy, 1);
    else
        % �������Ϊ D �����Ʋ��ԣ�A �����Ʋ���
        % ��ôѡ�� A ���Ʋ����£�D ������ߵĲ���
        A_strategy_index = find(dominant_A_strategy, 1);
        D_subopt_benefit_matrix = expected_D_matrix(:,A_strategy_index);
        [~,D_opt_strategy] = max(D_subopt_benefit_matrix);
        % �����Ʋ��Զ�Ӧ����������dominant_D_strategy
        dominant_D_strategy(D_opt_strategy) = 1;
    end
    if any(dominant_A_strategy)
        A_opt_strategy = find(dominant_A_strategy, 1);
    else
        % �������Ϊ A �����Ʋ��ԣ�D �����Ʋ���
        % ��ôѡ�� D ���Ʋ����£�A ������ߵĲ���
        D_strategy_index = find(dominant_D_strategy, 1);
        A_subopt_benefit_matrix = expected_A_matrix(D_strategy_index,:);
        [~,A_opt_strategy] = max(A_subopt_benefit_matrix);
        % �����Ʋ��Զ�Ӧ����������dominant_A_strategy
        dominant_A_strategy(A_opt_strategy) = 1;
    end

    % ����������������
    efficiency_D = zeros(1,nD);
    efficiency_A = zeros(1,nA);
    F_star_A = dominant_A_strategy;
    F_star_D = dominant_D_strategy;
    for k = 1:nD
        efficiency_D(k) = calculate_efficiency(num_type, P_A, F_star_A, U, k, 'defender');
    end
    for k = 1:nA
        efficiency_A(k) = calculate_efficiency(num_type, P_D, F_star_D, U, k, 'attacker');
    end

    % ��������ȷ���������Դ�������
    A_strategies = cell(1, nA);
    D_strategies = cell(1, nD);
    for i = 1:nA
        A_strategies{i} = sprintf('a%d', i);
    end
    for i = 1:nD
        D_strategies{i} = sprintf('d%d', i);
    end
    optimal_A_strategy = A_strategies{A_opt_strategy};
    optimal_D_strategy = D_strategies{D_opt_strategy};

    disp(['A �������Ų���: ', num2str(optimal_A_strategy)]);
    disp(['A ���Ĺ�������Ч��������',num2str(efficiency_A)]);
    disp(['D �������Ų���: ', num2str(optimal_D_strategy)]);
    disp(['D ���Ĺ�������Ч��������',num2str(efficiency_D)]);
end

%% �Ӻ�������
function efficiency = calculate_efficiency_ver2(num_type, P, F_star, U, strategy, player)
% ˵����������������Ч��
% ���㹫ʽ������̬��Ҷ˹����������������ѡȡ������ʽ(1)
% ����˵����
% ���룺
% num_type ���� [������, ������] ������
% P ���� �ҷ��Եз����͵��������
% F_star ���� �з�����ѡȡ����
% U ���� �����������
% strategy ���� �ҷ�ѡȡ�ĵ�i������
% player ���� �ҷ���ݣ�������defender���͡�attacker��
% �����
% efficiency ���� ��������Ч��
% �汾��2.0 ���� ��F_star(m_A_i)�޸�ΪF_star(t_A,m_A_i) ��2025/02/25��
%--------------------------------------------------------------------------
    efficiency = 0;
    num_type_A = num_type(1);   % ������������
    num_type_D = num_type(2);   % ������������
    if strcmp(player,'defender')
        nA = length(F_star);    % Attacker�Ĳ�����
        player_mode = 1;        % �� U ��ʹ�ã�ȡ�������ĵ� 1 λ��
        for t_A = 1:num_type_A  % �������п��ܵĹ��������� t_A
            for m_A_i = 1:nA    % �������п��ܵĹ����߲��� m_A_i
                % ���㲢�ۼ�Ч������ֵ������ÿ�ֹ��������ͺ��ж�����ϣ���������� U �л�ȡ��Ӧ��Ч�ú���ֵ��
                % ���Թ��������͵ĸ��� P((t_A, t_D) �͹��������͵ĸ��� F_star_A(t_A)��Ȼ�󽫽���ӵ� efficiency �����ϡ�
                efficiency = efficiency + P(t_A) * U{strategy, m_A_i, t_A}(player_mode) * F_star(t_A,m_A_i);
            end
        end                
    elseif strcmp(player,'attacker')
        nD = length(F_star);    % Defender�Ĳ�����      
        player_mode = 2;        % �� U ��ʹ�ã�ȡ�������ĵ� 2 λ��
        for t_A = 1:num_type_A  % �������п��ܵĹ��������� t_A
            for s_D_j = 1:nD    % �������п��ܵķ����߲��� m_A_i
                % **** ����Ҫ�� **** ��ΪF_starҪ�÷����ߵģ�״̬��Done��ʱ�䣺2025/02/20��
                % ���㲢�ۼ�Ч������ֵ������ÿ�ֹ��������ͺ��ж�����ϣ���������� U �л�ȡ��Ӧ��Ч�ú���ֵ��
                % ���Թ��������͵ĸ��� P((t_A, t_D) �͹��������͵ĸ��� F_star_A(t_A)��Ȼ�󽫽���ӵ� efficiency �����ϡ�
                efficiency = efficiency + P(num_type_D) * U{s_D_j, strategy, t_A}(player_mode) * F_star(s_D_j);
            end
        end
    end
end