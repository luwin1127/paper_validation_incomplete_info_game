clc;clear;close all;

% 论文《Near-optimal evalutation of network survivability under multi-stage attacks》扩展为不完全信息静态博弈
% 定义概率，攻击者是强类型的概率
p = 0.6;

% 定义强类型攻击者下的收益矩阵，使用元胞数组存储
DOD_matrix_strong = {
    [2, 8], [1, 9], [0.5, 9.5], [0.2, 9.8], [0.1, 9.9];
    [2.5, 7.5], [2, 8], [1.5, 8.5], [1.2, 8.8], [1.1, 8.9];
    [3, 7], [2.5, 7.5], [2, 8], [1.8, 8.2], [1.7, 8.3];
    [3.5, 6.5], [3, 7], [2.5, 7.5], [2.2, 7.8], [2.1, 7.9];
    [4, 6], [3.5, 6.5], [3, 7], [2.8, 7.2], [2.7, 7.3]
};

% 定义弱类型攻击者下的收益矩阵，使用元胞数组存储
DOD_matrix_weak = {
    [3, 7], [2, 8], [1.5, 8.5], [1.2, 8.8], [1.1, 8.9];
    [3.5, 6.5], [3, 7], [2.5, 7.5], [2.2, 7.8], [2.1, 7.9];
    [4, 6], [3.5, 6.5], [3, 7], [2.8, 7.2], [2.7, 7.3];
    [4.5, 5.5], [4, 6], [3.5, 6.5], [3.2, 6.8], [3.1, 6.9];
    [5, 5], [4.5, 5.5], [4, 6], [3.8, 6.2], [3.7, 6.3]
};

% 获取防御者策略的数量（矩阵的行数）
num_defender_strategies = size(DOD_matrix_strong, 1);
% 获取攻击者策略的数量（矩阵的列数）
num_attacker_strategies = size(DOD_matrix_strong, 2);

% 初始化期望收益矩阵，分别存储防御者和攻击者的期望收益
expected_defender_matrix = zeros(num_defender_strategies, num_attacker_strategies);
expected_attacker_matrix = zeros(num_defender_strategies, num_attacker_strategies);

% 计算期望收益矩阵
for i = 1:num_defender_strategies
    for j = 1:num_attacker_strategies
        expected_defender_matrix(i, j) = p * DOD_matrix_strong{i, j}(1) + (1 - p) * DOD_matrix_weak{i, j}(1);
        expected_attacker_matrix(i, j) = p * DOD_matrix_strong{i, j}(2) + (1 - p) * DOD_matrix_weak{i, j}(2);
    end
end

% 检查优势策略部分 - 防御者
% 初始化一个布尔向量，用于标记防御者的每个策略是否为优势策略，初始都设为 false
dominant_defender_strategy = false(num_defender_strategies, 1);
% 外层循环，遍历防御者的每一个策略
for i = 1:num_defender_strategies
    % 先假设当前策略 i 是优势策略
    is_dominant = true;
    % 中间层循环，将当前策略 i 与其他所有防御者策略进行比较
    for j = 1:num_defender_strategies
        % 跳过与自身的比较
        if i ~= j
            % 内层循环，对于攻击者的每一个策略
            for k = 1:num_attacker_strategies
                % 如果在某个攻击者策略 k 下，当前策略 i 的期望收益小于策略 j 的期望收益
                % 说明当前策略 i 不是优势策略
                if expected_defender_matrix(i, k) < expected_defender_matrix(j, k)
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
    dominant_defender_strategy(i) = is_dominant;
end

% 检查优势策略部分 - 攻击者
% 初始化一个布尔向量，用于标记攻击者的每个策略是否为优势策略，初始都设为 false
dominant_attacker_strategy = false(num_attacker_strategies, 1);
% 外层循环，遍历攻击者的每一个策略
for i = 1:num_attacker_strategies
    % 先假设当前策略 i 是优势策略
    is_dominant = true;
    % 中间层循环，将当前策略 i 与其他所有攻击者策略进行比较
    for j = 1:num_attacker_strategies
        % 跳过与自身的比较
        if i ~= j
            % 内层循环，对于防御者的每一个策略
            for k = 1:num_defender_strategies
                % 如果在某个防御者策略 k 下，当前策略 i 的期望收益小于策略 j 的期望收益
                % 说明当前策略 i 不是优势策略
                if expected_attacker_matrix(k, i) < expected_attacker_matrix(k, j)
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
    dominant_attacker_strategy(i) = is_dominant;
end

% 如果没有优势策略，使用期望收益矩阵进行求解
if ~any(dominant_defender_strategy) && ~any(dominant_attacker_strategy)
    % 初始化攻击者和防御者的最优策略索引
    attacker_opt_strategy = 1;
    defender_opt_strategy = 1;

    % 初始化最大期望收益
    max_attacker_income = -Inf;
    max_defender_income = -Inf;

    % 寻找攻击者的最优策略（最大化期望收益）
    for i = 1:num_attacker_strategies
        total_income = sum(expected_attacker_matrix(:, i));
        if total_income > max_attacker_income
            max_attacker_income = total_income;
            attacker_opt_strategy = i;
        end
    end

    % 寻找防御者的最优策略（最大化期望收益）
    for j = 1:num_defender_strategies
        total_income = sum(expected_defender_matrix(j, :));
        if total_income > max_defender_income
            max_defender_income = total_income;
            defender_opt_strategy = j;
        end
    end
else
    % 如果有优势策略，直接根据优势策略确定最优策略
    if any(dominant_defender_strategy)
        defender_opt_strategy = find(dominant_defender_strategy, 1);
    end
    if any(dominant_attacker_strategy)
        attacker_opt_strategy = find(dominant_attacker_strategy, 1);
    end
end

% 根据索引确定具体的资源分配策略
attacker_strategies = {[0, 20], [5, 15], [10, 10], [15, 5], [20, 0]};
defender_strategies = {[0, 20], [5, 15], [10, 10], [15, 5], [20, 0]};
optimal_attacker_strategy = attacker_strategies{attacker_opt_strategy};
optimal_defender_strategy = defender_strategies{defender_opt_strategy};

disp(['攻击者的最优策略: ', num2str(optimal_attacker_strategy)]);
disp(['防御者的最优策略: ', num2str(optimal_defender_strategy)]);