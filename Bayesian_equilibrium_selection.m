clc;clear;close all;

% 论文《静态贝叶斯博弈主动防御策略选取方法》验证
% 定义先验概率
P_D = [1,0];        % 攻击者对【防御者】类型的先验信念
P_A = [0.5,0.5];    % 防御者对【攻击者】类型的先验信念

% 混合策略
F_star_D = [0.5, 0.5];
F_star_A = [0.2, 0.8];

% 定义收益矩阵，每个元素 (m_A_i, s_D_j, t_A) 对应一个收益值
% 每个博弈收益对中，第 1 个为防御者的博弈收益，第 2 个数为攻击者的博弈收益
U = cat(3, ...
    {[40,50], [-10,0]; [0,300], [0,300]}, ...   % 攻击者类型1
    {[30,80], [-10,100]; [0,400], [0,400]});    % 攻击者类型2

% 防御策略效能
efficiency_D1 = calculate_efficiency(P_D, F_star_D, U, 1, 'defender');
efficiency_D2 = calculate_efficiency(P_D, F_star_D, U, 2, 'defender');

fprintf('E_s_D1: %.4f\n', efficiency_D1);
fprintf('E_s_D2: %.4f\n', efficiency_D2);

% 攻击策略效能
efficiency_A1 = calculate_efficiency(P_A, F_star_A, U, 1, 'attacker');
efficiency_A2 = calculate_efficiency(P_A, F_star_A, U, 2, 'attacker');

fprintf('E_s_A1: %.4f\n', efficiency_A1);
fprintf('E_s_A2: %.4f\n', efficiency_A2);

function efficiency = calculate_efficiency(P, F_star, U, strategy, player)
% 说明：计算期望博弈效能
% 计算公式：《静态贝叶斯博弈主动防御策略选取方法》式(1)
% 变量说明：
% 输入：
% P —— 我方对敌方类型的先验概率
% F_star —— 敌方策略选取概率
% U —— 博弈收益矩阵
% strategt —— 我方选取的第i个策略
% player —— 我方身份，包括“defender”和“attacker”
% 输出：
% efficiency —— 期望博弈效能
%--------------------------------------------------------------------------
    efficiency = 0;
    if strcmp(player,'defender')
        nD = length(strategy);  % Defender的策略数
        nA = length(F_star);    % Attacker的策略数
        player_mode = 1;        % 在 U 里使用，取收益矩阵的第 1 位数
        for t_A = 1:nA          % 遍历所有可能的攻击者类型 t_A
            for m_A_i = 1:nD    % 遍历所有可能的攻击者策略 m_A_i
                % 计算并累加效能期望值。对于每种攻击者类型和行动的组合，从收益矩阵 U 中获取相应的效用函数值，
                % 乘以攻击者类型的概率 P((t_A, t_D) 和攻击者类型的概率 F_star_A(t_A)，然后将结果加到 efficiency 变量上。
                efficiency = efficiency + P(t_A) * U{m_A_i, strategy, t_A}(player_mode) * F_star(m_A_i);
            end
        end                
    elseif strcmp(player,'attacker')
        nD = length(strategy);  % Defender的策略数
        nA = length(F_star);    % Attacker的策略数        
        player_mode = 2;        % 在 U 里使用，取收益矩阵的第 2 位数
        for t_A = 1:nA          % 遍历所有可能的攻击者类型 t_A
            for s_D_j = 1:nD    % 遍历所有可能的防守者策略 m_A_i
                % **** 这里要改 **** 因为F_star要用防御者的（状态：Done；时间：2025/02/20）
                % 计算并累加效能期望值。对于每种攻击者类型和行动的组合，从收益矩阵 U 中获取相应的效用函数值，
                % 乘以攻击者类型的概率 P((t_A, t_D) 和攻击者类型的概率 F_star_A(t_A)，然后将结果加到 efficiency 变量上。
                efficiency = efficiency + P(1) * U{strategy, s_D_j, t_A}(player_mode) * F_star(s_D_j);
            end
        end
    end
end
