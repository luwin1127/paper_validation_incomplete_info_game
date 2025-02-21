clc;clear;close all;

% 论文《Near-optimal evalutation of network survivability under multi-stage attacks》验证
% 定义收益矩阵，这里根据论文中的数据构建
% 行表示防御者策略，列表示攻击者策略
DOD_matrix = [
    0.749355, 2.979505, 2.99328, 2.995515, 2.996635;
    0.854925, 1.870855, 2.305035, 2.448675, 2.449975;
    0.99734, 1.632309, 1.872685, 2.0907, 2.148695;
    1.197075, 1.624254, 1.821465, 1.872965, 1.84314;
    1.496635, 1.798021, 1.997645, 2.13829, 1.49871
];

% 获取防御者策略的数量（矩阵的行数）
num_defender_strategies = size(DOD_matrix, 1);
% 获取攻击者策略的数量（矩阵的列数）
num_attacker_strategies = size(DOD_matrix, 2);

% 用于标记防御者的每个策略是否为优势策略，初始都设为 false
dominant_defender_strategy = false(num_defender_strategies, 1);
% 遍历防御者的每一个策略
for i = 1:num_defender_strategies
    % 先假设当前策略 i 是优势策略
    is_dominant = true;
    % 将当前策略 i 与其他所有防御者策略进行比较
    for j = 1:num_defender_strategies
        % 跳过与自身的比较
        if i ~= j
            % 对于攻击者的每一个策略
            for k = 1:num_attacker_strategies
                % 如果在某个攻击者策略 k 下，当前策略 i 的 DOD 值大于策略 j 的 DOD 值
                % 说明当前策略 i 不是优势策略
                if DOD_matrix(i, k) > DOD_matrix(j, k)
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

% 用于标记攻击者的每个策略是否为优势策略，初始都设为 false
dominant_attacker_strategy = false(num_attacker_strategies, 1);
% 遍历攻击者的每一个策略
for i = 1:num_attacker_strategies
    % 先假设当前策略 i 是优势策略
    is_dominant = true;
    % 将当前策略 i 与其他所有攻击者策略进行比较
    for j = 1:num_attacker_strategies
        % 跳过与自身的比较
        if i ~= j
            % 对于防御者的每一个策略
            for k = 1:num_defender_strategies
                % 如果在某个防御者策略 k 下，当前策略 i 的 DOD 值小于策略 j 的 DOD 值
                % 说明当前策略 i 不是优势策略
                if DOD_matrix(k, i) < DOD_matrix(k, j)
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

% 如果没有优势策略，使用原矩阵进行最小最大策略求解
if ~any(dominant_defender_strategy) && ~any(dominant_attacker_strategy)
    % 初始化攻击者和防御者的最优策略索引
    attacker_opt_strategy = 1;
    defender_opt_strategy = 1;

    % 初始化最小最大和最大最小值
    max_min_attacker = -Inf;
    min_max_defender = Inf;

    % 寻找攻击者的最优策略（最大化最小收益）
    for i = 1:num_attacker_strategies
        min_value = min(DOD_matrix(:, i));
        if min_value > max_min_attacker
            max_min_attacker = min_value;
            attacker_opt_strategy = i;
        end
    end

    % 寻找防御者的最优策略（最小化最大损失）
    for j = 1:num_defender_strategies
        max_value = max(DOD_matrix(j, :));
        if max_value < min_max_defender
            min_max_defender = max_value;
            defender_opt_strategy = j;
        end
    end
end

% 根据索引确定具体的资源分配策略
attacker_strategies = {[0, 20], [5, 15], [10, 10], [15, 5], [20, 0]};
defender_strategies = {[0, 20], [5, 15], [10, 10], [15, 5], [20, 0]};
optimal_attacker_strategy = attacker_strategies{attacker_opt_strategy};
optimal_defender_strategy = defender_strategies{defender_opt_strategy};

disp(['攻击者的最优策略: ', num2str(optimal_attacker_strategy)]);
disp(['防御者的最优策略: ', num2str(optimal_defender_strategy)]);