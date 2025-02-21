clc;clear;close all;

% 书籍《博弈论与信息经济学》pp. 143-144 市场阻挠博弈：不完全信息
% 说明：读完《Near-optimal evaluation of network survivability under multi-stage
% attacks》和《Resource allocation strategies to maximize network survivability 
% considering of average DOD》
% 求解过程分为2步：
% 1. 构建收益矩阵
% 2. 寻找优势策略
% （1）优势策略消除
% （2）最小最大策略
% （3）混合策略 线性规划

% 定义概率，在位者是高成本的概率
p = 0.3;

% 定义均衡求解方式
mode = 'mixed';     % 混合策略
% mode = 'pure';      % 纯策略

% 定义高成本情况在位者下的收益矩阵，使用元胞数组存储
high_benefit_matrix = {
  [40, 50], [-10, 0];
  [0, 300], [0, 300]
};

% 定义低成本情况在位者下的收益矩阵，使用元胞数组存储
low_benefit_matrix = {
  [30, 80], [-10, 100];
  [0, 400], [0, 400]
};

% 获取进入者（entrant）策略的数量（矩阵的行数）
num_entrant_strategies = size(high_benefit_matrix, 1);
% 获取在位者（incumbent）策略的数量（矩阵的列数）
num_incumbent_strategies = size(high_benefit_matrix, 2);

% 初始化期望收益矩阵，分别存储进入者和在位者的期望收益
expected_entrant_matrix = zeros(num_entrant_strategies, num_incumbent_strategies);
expected_incumbent_matrix = zeros(num_entrant_strategies, num_incumbent_strategies);

% 计算期望收益矩阵
for i = 1:num_entrant_strategies
    for j = 1:num_incumbent_strategies
        expected_entrant_matrix(i, j) = p * high_benefit_matrix{i, j}(1) + (1 - p) * low_benefit_matrix{i, j}(1);
        expected_incumbent_matrix(i, j) = p * high_benefit_matrix{i, j}(2) + (1 - p) * low_benefit_matrix{i, j}(2);
    end
end

% 检查优势策略部分 - 进入者
% 初始化一个布尔向量，用于标记进入者的每个策略是否为优势策略，初始都设为 false
dominant_entrant_strategy = false(num_entrant_strategies, 1);
% 外层循环，遍历进入者的每一个策略
for i = 1:num_entrant_strategies
    % 先假设当前策略 i 是优势策略
    is_dominant = true;
    % 中间层循环，将当前策略 i 与其他所有进入者策略进行比较
    for j = 1:num_entrant_strategies
        % 跳过与自身的比较
        if i ~= j
            % 内层循环，对于在位者的每一个策略
            for k = 1:num_incumbent_strategies
                % 如果在某个在位者策略 k 下，当前策略 i 的期望收益小于策略 j 的期望收益
                % 说明当前策略 i 不是优势策略
                if expected_entrant_matrix(i, k) < expected_entrant_matrix(j, k)
                    is_dominant = false;
                    % 一旦发现不是优势策略，就跳出内层循环
                    break;
                end
            end
            % 如果已经确定不是优势策略，跳出中间层循环
            if ~is_dominant
                break;
            end
        end
    end
    % 如果经过所有比较后，当前策略 i 仍然被认为是优势策略
    % 则将对应的标记位置为 true
    dominant_entrant_strategy(i) = is_dominant;
end

% 检查优势策略部分 - 在位者
% 初始化一个布尔向量，用于标记在位者的每个策略是否为优势策略，初始都设为 false
dominant_incumbent_strategy = false(num_incumbent_strategies, 1);
% 外层循环，遍历在位者的每一个策略
for i = 1:num_incumbent_strategies
    % 先假设当前策略 i 是优势策略
    is_dominant = true;
    % 中间层循环，将当前策略 i 与其他所有在位者策略进行比较
    for j = 1:num_incumbent_strategies
        % 跳过与自身的比较
        if i ~= j
            % 内层循环，对于进入者的每一个策略
            for k = 1:num_entrant_strategies
                % 如果在某个进入者策略 k 下，当前策略 i 的期望收益小于策略 j 的期望收益
                % 说明当前策略 i 不是优势策略
                if expected_incumbent_matrix(k, i) < expected_incumbent_matrix(k, j)
                    is_dominant = false;
                    % 一旦发现不是优势策略，就跳出内层循环
                    break;
                end
            end
            % 如果已经确定不是优势策略，跳出中间层循环
            if ~is_dominant
                break;
            end
        end
    end
    % 如果经过所有比较后，当前策略 i 仍然被认为是优势策略
    % 则将对应的标记位置为 true
    dominant_incumbent_strategy(i) = is_dominant;
end

% 如果没有优势策略，使用期望收益矩阵进行求解
if ~any(dominant_entrant_strategy) && ~any(dominant_incumbent_strategy)
    if strcmp(mode,'mixed')
        % 进入者混合策略求解
        % 目标函数：最小化一个虚拟变量 v（表示最小期望收益）
        f_entrant = [1, zeros(1, num_entrant_strategies)]; 
        % 约束条件：期望收益大于等于 v
        A_entrant = [-ones(num_incumbent_strategies, 1), expected_entrant_matrix];
        b_entrant = zeros(num_incumbent_strategies, 1);
        % 概率和为 1
        Aeq_entrant = [0, ones(1, num_entrant_strategies)];
        beq_entrant = 1;
        % 概率非负
        lb_entrant = [-Inf, zeros(1, num_entrant_strategies)];
        ub_entrant = [Inf, ones(1, num_entrant_strategies)];
        % 求解线性规划
        options = optimoptions('linprog', 'Display', 'off');
        [x_entrant, ~] = linprog(f_entrant, A_entrant, b_entrant, Aeq_entrant, beq_entrant, lb_entrant, ub_entrant, options);
        entrant_mixed_strategy = x_entrant(2:end);
        entrant_min_expected_payoff = x_entrant(1);

        % 在位者混合策略求解
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

        disp('进入者的混合策略:');
        disp(entrant_mixed_strategy);
        disp(['进入者的最小期望收益: ', num2str(entrant_min_expected_payoff)]);
        disp('在位者的混合策略:');
        disp(incumbent_mixed_strategy);
        disp(['在位者的最小期望收益: ', num2str(incumbent_min_expected_payoff)]);
    elseif strcmp(mode,'pure')
        % 初始化在位者和进入者的最优策略索引
        incumbent_opt_strategy = 1;
        entrant_opt_strategy = 1;

        % 初始化最大期望收益
        max_incumbent_income = -Inf;
        max_entrant_income = -Inf;

        % 寻找在位者的最优策略（最大化期望收益）
        for i = 1:num_incumbent_strategies
            total_income = sum(expected_incumbent_matrix(:, i));
            if total_income > max_incumbent_income
                max_incumbent_income = total_income;
                incumbent_opt_strategy = i;
            end
        end

        % 寻找进入者的最优策略（最大化期望收益）
        for j = 1:num_entrant_strategies
            total_income = sum(expected_entrant_matrix(j, :));
            if total_income > max_entrant_income
                max_entrant_income = total_income;
                entrant_opt_strategy = j;
            end
        end
    end
else
    % 如果有优势策略，直接根据优势策略确定最优策略
    if any(dominant_entrant_strategy)
        entrant_opt_strategy = find(dominant_entrant_strategy, 1);
    else
        % 这种情况为进入者无优势策略，在位者有优势策略
        % 那么选择在位者优势策略下，进入者收益最高的策略
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
        % 这种情况为在位者无优势策略，进入者有优势策略
        % 那么选择进入者优势策略下，在位者收益最高的策略
        entrant_strategy_index = find(dominant_incumbent_strategy, 1);
        if expected_entrant_matrix(1,entrant_strategy_index) >= expected_entrant_matrix(2, entrant_strategy_index)
            incumbent_opt_strategy = 1;
        else
            incumbent_opt_strategy = 2;
        end
    end
end

% 根据索引确定具体的资源分配策略
incumbent_strategies = {"默许","斗争"};
entrant_strategies = {"进入","不进入"};
optimal_incumbent_strategy = incumbent_strategies{incumbent_opt_strategy};
optimal_entrant_strategy = entrant_strategies{entrant_opt_strategy};

disp(['在位者的最优策略: ', num2str(optimal_incumbent_strategy)]);
disp(['进入者的最优策略: ', num2str(optimal_entrant_strategy)]);