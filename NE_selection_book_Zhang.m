clc;clear;close all;

% �鼮������������Ϣ����ѧ��pp. 143-144 �г����Ӳ��ģ�����ȫ��Ϣ
% ˵�������꡶Near-optimal evaluation of network survivability under multi-stage
% attacks���͡�Resource allocation strategies to maximize network survivability 
% considering of average DOD��
% �����̷�Ϊ2����
% 1. �����������
% 2. Ѱ�����Ʋ���
% ��1�����Ʋ�������
% ��2����С������
% ��3����ϲ��� ���Թ滮

% ������ʣ���λ���Ǹ߳ɱ��ĸ���
p = 0.3;

% ���������ⷽʽ
mode = 'mixed';     % ��ϲ���
% mode = 'pure';      % ������

% ����߳ɱ������λ���µ��������ʹ��Ԫ������洢
high_benefit_matrix = {
  [40, 50], [-10, 0];
  [0, 300], [0, 300]
};

% ����ͳɱ������λ���µ��������ʹ��Ԫ������洢
low_benefit_matrix = {
  [30, 80], [-10, 100];
  [0, 400], [0, 400]
};

% ��ȡ�����ߣ�entrant�����Ե������������������
num_entrant_strategies = size(high_benefit_matrix, 1);
% ��ȡ��λ�ߣ�incumbent�����Ե������������������
num_incumbent_strategies = size(high_benefit_matrix, 2);

% ��ʼ������������󣬷ֱ�洢�����ߺ���λ�ߵ���������
expected_entrant_matrix = zeros(num_entrant_strategies, num_incumbent_strategies);
expected_incumbent_matrix = zeros(num_entrant_strategies, num_incumbent_strategies);

% ���������������
for i = 1:num_entrant_strategies
    for j = 1:num_incumbent_strategies
        expected_entrant_matrix(i, j) = p * high_benefit_matrix{i, j}(1) + (1 - p) * low_benefit_matrix{i, j}(1);
        expected_incumbent_matrix(i, j) = p * high_benefit_matrix{i, j}(2) + (1 - p) * low_benefit_matrix{i, j}(2);
    end
end

% ������Ʋ��Բ��� - ������
% ��ʼ��һ���������������ڱ�ǽ����ߵ�ÿ�������Ƿ�Ϊ���Ʋ��ԣ���ʼ����Ϊ false
dominant_entrant_strategy = false(num_entrant_strategies, 1);
% ���ѭ�������������ߵ�ÿһ������
for i = 1:num_entrant_strategies
    % �ȼ��赱ǰ���� i �����Ʋ���
    is_dominant = true;
    % �м��ѭ��������ǰ���� i ���������н����߲��Խ��бȽ�
    for j = 1:num_entrant_strategies
        % ����������ıȽ�
        if i ~= j
            % �ڲ�ѭ����������λ�ߵ�ÿһ������
            for k = 1:num_incumbent_strategies
                % �����ĳ����λ�߲��� k �£���ǰ���� i ����������С�ڲ��� j ����������
                % ˵����ǰ���� i �������Ʋ���
                if expected_entrant_matrix(i, k) < expected_entrant_matrix(j, k)
                    is_dominant = false;
                    % һ�����ֲ������Ʋ��ԣ��������ڲ�ѭ��
                    break;
                end
            end
            % ����Ѿ�ȷ���������Ʋ��ԣ������м��ѭ��
            if ~is_dominant
                break;
            end
        end
    end
    % ����������бȽϺ󣬵�ǰ���� i ��Ȼ����Ϊ�����Ʋ���
    % �򽫶�Ӧ�ı��λ��Ϊ true
    dominant_entrant_strategy(i) = is_dominant;
end

% ������Ʋ��Բ��� - ��λ��
% ��ʼ��һ���������������ڱ����λ�ߵ�ÿ�������Ƿ�Ϊ���Ʋ��ԣ���ʼ����Ϊ false
dominant_incumbent_strategy = false(num_incumbent_strategies, 1);
% ���ѭ����������λ�ߵ�ÿһ������
for i = 1:num_incumbent_strategies
    % �ȼ��赱ǰ���� i �����Ʋ���
    is_dominant = true;
    % �м��ѭ��������ǰ���� i ������������λ�߲��Խ��бȽ�
    for j = 1:num_incumbent_strategies
        % ����������ıȽ�
        if i ~= j
            % �ڲ�ѭ�������ڽ����ߵ�ÿһ������
            for k = 1:num_entrant_strategies
                % �����ĳ�������߲��� k �£���ǰ���� i ����������С�ڲ��� j ����������
                % ˵����ǰ���� i �������Ʋ���
                if expected_incumbent_matrix(k, i) < expected_incumbent_matrix(k, j)
                    is_dominant = false;
                    % һ�����ֲ������Ʋ��ԣ��������ڲ�ѭ��
                    break;
                end
            end
            % ����Ѿ�ȷ���������Ʋ��ԣ������м��ѭ��
            if ~is_dominant
                break;
            end
        end
    end
    % ����������бȽϺ󣬵�ǰ���� i ��Ȼ����Ϊ�����Ʋ���
    % �򽫶�Ӧ�ı��λ��Ϊ true
    dominant_incumbent_strategy(i) = is_dominant;
end

% ���û�����Ʋ��ԣ�ʹ�������������������
if ~any(dominant_entrant_strategy) && ~any(dominant_incumbent_strategy)
    if strcmp(mode,'mixed')
        % �����߻�ϲ������
        % Ŀ�꺯������С��һ��������� v����ʾ��С�������棩
        f_entrant = [1, zeros(1, num_entrant_strategies)]; 
        % Լ������������������ڵ��� v
        A_entrant = [-ones(num_incumbent_strategies, 1), expected_entrant_matrix];
        b_entrant = zeros(num_incumbent_strategies, 1);
        % ���ʺ�Ϊ 1
        Aeq_entrant = [0, ones(1, num_entrant_strategies)];
        beq_entrant = 1;
        % ���ʷǸ�
        lb_entrant = [-Inf, zeros(1, num_entrant_strategies)];
        ub_entrant = [Inf, ones(1, num_entrant_strategies)];
        % ������Թ滮
        options = optimoptions('linprog', 'Display', 'off');
        [x_entrant, ~] = linprog(f_entrant, A_entrant, b_entrant, Aeq_entrant, beq_entrant, lb_entrant, ub_entrant, options);
        entrant_mixed_strategy = x_entrant(2:end);
        entrant_min_expected_payoff = x_entrant(1);

        % ��λ�߻�ϲ������
        f_incumbent = [1, zeros(1, num_incumbent_strategies)];
        A_incumbent = [-ones(num_entrant_strategies, 1), expected_incumbent_matrix'];
        b_incumbent = zeros(num_entrant_strategies, 1);
        Aeq_incumbent = [0, ones(1, num_incumbent_strategies)];
        beq_incumbent = 1;
        lb_incumbent = [-Inf, zeros(1, num_incumbent_strategies)];
        ub_incumbent = [Inf, ones(1, num_incumbent_strategies)];
        [x_incumbent, ~] = linprog(f_incumbent, A_incumbent, b_incumbent, Aeq_incumbent, beq_incumbent, lb_incumbent, ub_incumbent, options);
        incumbent_mixed_strategy = x_incumbent(2:end);
        incumbent_min_expected_payoff = x_incumbent(1);

        disp('�����ߵĻ�ϲ���:');
        disp(entrant_mixed_strategy);
        disp(['�����ߵ���С��������: ', num2str(entrant_min_expected_payoff)]);
        disp('��λ�ߵĻ�ϲ���:');
        disp(incumbent_mixed_strategy);
        disp(['��λ�ߵ���С��������: ', num2str(incumbent_min_expected_payoff)]);
    elseif strcmp(mode,'pure')
        % ��ʼ����λ�ߺͽ����ߵ����Ų�������
        incumbent_opt_strategy = 1;
        entrant_opt_strategy = 1;

        % ��ʼ�������������
        max_incumbent_income = -Inf;
        max_entrant_income = -Inf;

        % Ѱ����λ�ߵ����Ų��ԣ�����������棩
        for i = 1:num_incumbent_strategies
            total_income = sum(expected_incumbent_matrix(:, i));
            if total_income > max_incumbent_income
                max_incumbent_income = total_income;
                incumbent_opt_strategy = i;
            end
        end

        % Ѱ�ҽ����ߵ����Ų��ԣ�����������棩
        for j = 1:num_entrant_strategies
            total_income = sum(expected_entrant_matrix(j, :));
            if total_income > max_entrant_income
                max_entrant_income = total_income;
                entrant_opt_strategy = j;
            end
        end
    end
else
    % ��������Ʋ��ԣ�ֱ�Ӹ������Ʋ���ȷ�����Ų���
    if any(dominant_entrant_strategy)
        entrant_opt_strategy = find(dominant_entrant_strategy, 1);
    else
        % �������Ϊ�����������Ʋ��ԣ���λ�������Ʋ���
        % ��ôѡ����λ�����Ʋ����£�������������ߵĲ���
        incumbent_strategy_index = find(dominant_incumbent_strategy, 1);
        if expected_entrant_matrix(1, incumbent_strategy_index) >= expected_entrant_matrix(2, incumbent_strategy_index)
            entrant_opt_strategy = 1;
        else
            entrant_opt_strategy = 2;
        end
    end
    if any(dominant_incumbent_strategy)
        incumbent_opt_strategy = find(dominant_incumbent_strategy, 1);
    else
        % �������Ϊ��λ�������Ʋ��ԣ������������Ʋ���
        % ��ôѡ����������Ʋ����£���λ��������ߵĲ���
        entrant_strategy_index = find(dominant_incumbent_strategy, 1);
        if expected_entrant_matrix(1,entrant_strategy_index) >= expected_entrant_matrix(2, entrant_strategy_index)
            incumbent_opt_strategy = 1;
        else
            incumbent_opt_strategy = 2;
        end
    end
end

% ��������ȷ���������Դ�������
incumbent_strategies = {"Ĭ��","����"};
entrant_strategies = {"����","������"};
optimal_incumbent_strategy = incumbent_strategies{incumbent_opt_strategy};
optimal_entrant_strategy = entrant_strategies{entrant_opt_strategy};

disp(['��λ�ߵ����Ų���: ', num2str(optimal_incumbent_strategy)]);
disp(['�����ߵ����Ų���: ', num2str(optimal_entrant_strategy)]);