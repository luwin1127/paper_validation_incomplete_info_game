clc;clear;

% 验证论文《静态贝叶斯博弈主动防御策略选取方法》的仿真实验部分——linprog()方法
% 定义概率，攻击者是冒险型的概率
p = 0.6;
P_A = [p,1-p];      % 防御者对【攻击者】类型的先验信念
P_D = [0.3,0.4,0.4];        % 攻击者对【防御者】类型的先验信念
num_type_A = 2;     % 攻击者类型数
num_type_D = 1;     % 防御者类型数
num_type = [num_type_A, num_type_D];    % 第1个数是攻击者类型数，第2个数是防御者类型数

% 定义冒险型攻击者情况下的收益矩阵
high_benefit_matrix = {
    [-730,1050], [-680,980], [-580,950];
    [-750,1200], [-740,1000], [-680,1100];
    [-560,900], [-700,980], [-680,1000]
};

% 定义保守型攻击者情况下的收益矩阵
low_benefit_matrix = {
    [-780,1280],[-500,810],[-540,900];
    [-700,1100],[-670,1100],[-610,1080];
    [-510,1000],[-580,1020],[-580,1180];
};

% 定义收益矩阵，用于计算期望博弈收益
U = cat(3, high_benefit_matrix, low_benefit_matrix);

% 获取 D 的策略数量（矩阵的行数）
nD = size(high_benefit_matrix, 1);
% 获取 A 的策略数量（矩阵的列数）
nA = size(high_benefit_matrix, 2);

% 初始化策略概率
F_A_case1 = [0,0.3417,0.6583];
F_A_case2 = [0.8547,0,0.1453];
F_star_A = [F_A_case1;F_A_case2];
F_star_D = [0,0,1];

% 初始化期望收益矩阵，分别存储 D 和 A 的期望收益
expected_D_matrix = zeros(nD, nA);
expected_A_matrix = zeros(nD, nA);

% 计算期望收益矩阵
for t_A = 1:num_type_A  % 遍历所有可能的攻击者类型 t_A
    for i = 1:nD
        for j = 1:nA
            expected_D_matrix(i, j) = p * high_benefit_matrix{i, j}(1) * F_star_D(i) + (1 - p) * low_benefit_matrix{i, j}(1) * F_star_D(i);
            expected_A_matrix(i, j) = p * high_benefit_matrix{i, j}(2) * F_star_A(t_A,j) + (1 - p) * low_benefit_matrix{i, j}(2) * F_star_A(t_A,j);
        end
    end
end

% 计算期望博弈收益
efficiency_D = zeros(1,nD);
efficiency_A = zeros(1,nA);

% 更新均衡下的策略
for k = 1:nD % 防御者各个策略的博弈收益
    efficiency_D(k) = calculate_efficiency_ver2(num_type, P_A, F_star_A, U, k, 'defender');
end
for k = 1:nA % 攻击者各个策略的博弈收益
    efficiency_A(k) = calculate_efficiency_ver2(num_type, P_D, F_star_D, U, k, 'attacker');
end

% 输出结果（这个部分证明了王晋东论文的正确）
disp('攻击者策略的攻击效能:');
disp(efficiency_A);
disp('防御者策略的防御效能:');
disp(efficiency_D);

% 检查优势策略部分 - D
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

% 检查优势策略部分 - A
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

% 如果没有优势策略，使用期望收益矩阵进行求解
if ~any(dominant_D_strategy) && ~any(dominant_A_strategy)
    if strcmp(mode, 'mixed')
        % D 混合策略求解
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
            error('D 混合策略线性规划求解失败');
        end
        F_star_D = x_D(2:end);
        D_max_expected_payoff = -x_D(1);

        % A 混合策略求解
        f_A = [-1, zeros(1, nA)];
        A_A = [ones(nD, 1), expected_A_matrix'];
        b_A = zeros(nD, 1);
        Aeq_A = [0, ones(1, nA)];
        beq_A = 1;
        lb_A = [-inf, zeros(1, nA)];
        ub_A = [inf, ones(1, nA)];
        [x_A, ~, exitflag_A] = linprog(f_A, A_A, b_A, Aeq_A, beq_A, lb_A, ub_A, options);
        if exitflag_A ~= 1
            error('A 混合策略线性规划求解失败');
        end
        F_star_A = x_A(2:end);
        A_min_expected_payoff = -x_A(1);

        % 计算期望博弈收益
        efficiency_D = zeros(1,nD);
        efficiency_A = zeros(1,nA);
        for k = 1:nD
            efficiency_D(k) = calculate_efficiency_ver2(num_type, P_A, F_star_A, U, k, 'defender');
        end
        for k = 1:nA
            efficiency_A(k) = calculate_efficiency_ver2(num_type, P_D, F_star_D, U, k, 'attacker');
        end

        disp('D 的混合策略:');
        disp(F_star_D);
        disp(['D 的最大期望收益: ', num2str(D_max_expected_payoff)]);
        D_opt_strategy = find(F_star_D);
        disp('A 的混合策略:');
        disp(F_star_A);
        disp(['A 的最大期望收益: ', num2str(A_min_expected_payoff)]);
        A_opt_strategy = find(F_star_A);

        % 根据索引确定具体的资源分配策略
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

        disp(['A 方的最优策略: ', strjoin(optimal_A_strategy, ', '),'，混合策略概率：',num2str(F_star_A')]);
        disp(['A 方的攻击策略效能期望：',num2str(efficiency_A)]);
        disp(['D 方的最优策略: ', strjoin(optimal_D_strategy, ', '),'，混合策略概率：',num2str(F_star_D')]);
        disp(['D 方的攻击策略效能期望：',num2str(efficiency_D)]);
    elseif strcmp(mode, 'pure')
        % 初始化 A 和 D 的最优策略索引
        A_opt_strategy = 1;
        D_opt_strategy = 1;

        % 初始化最大期望收益
        max_A_income = -Inf;
        max_D_income = -Inf;

        % 寻找 A 的最优策略（最大化期望收益）
        for i = 1:nA
            total_income = sum(expected_A_matrix(:, i));
            if total_income > max_A_income
                max_A_income = total_income;
                A_opt_strategy = i;
            end
        end

        % 寻找 D 的最优策略（最大化期望收益）
        for j = 1:nD
            total_income = sum(expected_D_matrix(j, :));
            if total_income > max_D_income
                max_D_income = total_income;
                D_opt_strategy = j;
            end
        end

        % 计算期望博弈收益
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

        % 根据索引确定具体的资源分配策略
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

        disp(['A 方的最优策略: ', num2str(optimal_A_strategy)]);
        disp(['A 方的攻击策略效能期望：',num2str(efficiency_A)]);
        disp(['D 方的最优策略: ', num2str(optimal_D_strategy)]);
        disp(['D 方的攻击策略效能期望：',num2str(efficiency_D)]);
    end
else
    % 如果有优势策略，并直接根据优势策略确定最优策略
    if any(dominant_D_strategy)
        D_opt_strategy = find(dominant_D_strategy, 1);
    else
        % 这种情况为 D 无优势策略，A 有优势策略
        % 那么选择 A 优势策略下，D 收益最高的策略
        A_strategy_index = find(dominant_A_strategy, 1);
        D_subopt_benefit_matrix = expected_D_matrix(:,A_strategy_index);
        [~,D_opt_strategy] = max(D_subopt_benefit_matrix);
        % 把优势策略对应的索引赋给dominant_D_strategy
        dominant_D_strategy(D_opt_strategy) = 1;
    end
    if any(dominant_A_strategy)
        A_opt_strategy = find(dominant_A_strategy, 1);
    else
        % 这种情况为 A 无优势策略，D 有优势策略
        % 那么选择 D 优势策略下，A 收益最高的策略
        D_strategy_index = find(dominant_D_strategy, 1);
        A_subopt_benefit_matrix = expected_A_matrix(D_strategy_index,:);
        [~,A_opt_strategy] = max(A_subopt_benefit_matrix);
        % 把优势策略对应的索引赋给dominant_A_strategy
        dominant_A_strategy(A_opt_strategy) = 1;
    end

    % 计算期望博弈收益
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

    % 根据索引确定具体的资源分配策略
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

    disp(['A 方的最优策略: ', num2str(optimal_A_strategy)]);
    disp(['A 方的攻击策略效能期望：',num2str(efficiency_A)]);
    disp(['D 方的最优策略: ', num2str(optimal_D_strategy)]);
    disp(['D 方的攻击策略效能期望：',num2str(efficiency_D)]);
end

%% 子函数部分
function efficiency = calculate_efficiency_ver2(num_type, P, F_star, U, strategy, player)
% 说明：计算期望博弈效能
% 计算公式：《静态贝叶斯博弈主动防御策略选取方法》式(1)
% 变量说明：
% 输入：
% num_type —— [攻击者, 防御者] 类型数
% P —— 我方对敌方类型的先验概率
% F_star —— 敌方策略选取概率
% U —— 博弈收益矩阵
% strategy —— 我方选取的第i个策略
% player —— 我方身份，包括“defender”和“attacker”
% 输出：
% efficiency —— 期望博弈效能
% 版本：2.0 —— 将F_star(m_A_i)修改为F_star(t_A,m_A_i) （2025/02/25）
%--------------------------------------------------------------------------
    efficiency = 0;
    num_type_A = num_type(1);   % 攻击者类型数
    num_type_D = num_type(2);   % 防御者类型数
    if strcmp(player,'defender')
        nA = length(F_star);    % Attacker的策略数
        player_mode = 1;        % 在 U 里使用，取收益矩阵的第 1 位数
        for t_A = 1:num_type_A  % 遍历所有可能的攻击者类型 t_A
            for m_A_i = 1:nA    % 遍历所有可能的攻击者策略 m_A_i
                % 计算并累加效能期望值。对于每种攻击者类型和行动的组合，从收益矩阵 U 中获取相应的效用函数值，
                % 乘以攻击者类型的概率 P((t_A, t_D) 和攻击者类型的概率 F_star_A(t_A)，然后将结果加到 efficiency 变量上。
                efficiency = efficiency + P(t_A) * U{strategy, m_A_i, t_A}(player_mode) * F_star(t_A,m_A_i);
            end
        end                
    elseif strcmp(player,'attacker')
        nD = length(F_star);    % Defender的策略数      
        player_mode = 2;        % 在 U 里使用，取收益矩阵的第 2 位数
        for t_A = 1:num_type_A  % 遍历所有可能的攻击者类型 t_A
            for s_D_j = 1:nD    % 遍历所有可能的防守者策略 m_A_i
                % **** 这里要改 **** 因为F_star要用防御者的（状态：Done；时间：2025/02/20）
                % 计算并累加效能期望值。对于每种攻击者类型和行动的组合，从收益矩阵 U 中获取相应的效用函数值，
                % 乘以攻击者类型的概率 P((t_A, t_D) 和攻击者类型的概率 F_star_A(t_A)，然后将结果加到 efficiency 变量上。
                efficiency = efficiency + P(num_type_D) * U{s_D_j, strategy, t_A}(player_mode) * F_star(s_D_j);
            end
        end
    end
end