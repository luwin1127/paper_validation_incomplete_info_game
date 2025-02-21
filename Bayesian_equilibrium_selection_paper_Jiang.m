function Bayesian_equilibrium_selection_paper_Jiang
% 验证论文《基于攻防博弈模型的网络安全测评和最优主动防御》姜伟
clc;clear;close all;

%% 初始参数设置
% 论文《静态贝叶斯博弈主动防御策略选取方法》验证
% 定义先验概率
P_D = [1,0];        % 攻击者对【防御者】类型的先验信念
P_A = [0.9,0.1];    % 防御者对【攻击者】类型的先验信念
type_P_D = 1;       % 防御者类型数
type_P_A = 2;       % 攻击者类型数

% 策略数量
nD = 6; % 防御者策略数量
nA = 3; % 攻击者策略数量

% 混合策略
% F_star_D = ones(1,nD)*(1/nD);
% F_star_A = ones(1,nA)*(1/nA);
% 下面的结果是最优结果
F_star_D = [28/53, 0, 0, 0, 0, 25/36];
F_star_A = [67/159, 0, 92/159];

% 定义收益矩阵，每个元素 (s_A, s_D) 对应一个收益值
U = [810 1510 1670 2470 2650 2300
     810 2470 1670 2470 2470 2300
    2310 2310 2470 1670  970 2300];

% 防御策略效能
efficiency_D = calculate_efficiency_paper_Jiang(F_star_A, U, F_star_D, 'defender');
fprintf('E_s_D: %.4f\n', efficiency_D);

% 攻击策略效能
efficiency_A = calculate_efficiency_paper_Jiang(F_star_D, U, F_star_A, 'attacker');
fprintf('E_s_A1: %.4f\n', efficiency_A);

%% 算法开始
% 循环参数
% tol = 1;
options = optimset('Display','off',...  % off/iter/notify-detailed
                   'Algorithm','sqp',...% 算法选择
                   'TolCon',1e-6,...    % 约束违反度度容差（正标量），默认值为 1e-6
                   'TolX',1e-6,...      % 关于正标量 x 的终止容差，SQP算法的默认值为 1e-6
                   'TolFun',1e-6);      % 一阶最优性的终止容差（正标量），默认值为 1e-6             
lb = [0;0];
ub = [1;1];
Aeq = [1 1];
beq = 1;

% 策略选取算法
x_D = fmincon(@objective_D,F_star_D,[],[],Aeq,beq,lb,ub,[],options);
F_star_D = x_D;
disp(F_star_D);

x_A = fmincon(@objective_A,F_star_A,[],[],Aeq,beq,lb,ub,[],options);
F_star_A = x_A;
disp(F_star_A);

%% 结果显示
% 防御策略效能
efficiency_D1 = calculate_efficiency_paper_Jiang(F_star_A, U, 1, 'defender');
efficiency_D2 = calculate_efficiency_paper_Jiang(F_star_A, U, 2, 'defender');

fprintf('E_s_D1: %.4f\n', efficiency_D1);
fprintf('E_s_D2: %.4f\n', efficiency_D2);

% 攻击策略效能
efficiency_A1 = calculate_efficiency_paper_Jiang(F_star_D, U, 1, 'attacker');
efficiency_A2 = calculate_efficiency_paper_Jiang(F_star_D, U, 2, 'attacker');

fprintf('E_s_A1: %.4f\n', efficiency_A1);
fprintf('E_s_A2: %.4f\n', efficiency_A2);
%% 子函数部分
% 计算效能期望
function efficiency_res = calculate_efficiency_paper_Jiang(F_star, U, F_strategy, player)
    if strcmp(player,'defender')
        efficiency = zeros(1,nD);
        for s_D_j = 1:nD
            for s_A_i = 1:nA % 遍历所有可能的攻击者策略 m_A_i
                % 计算并累加效能期望值。对于每种攻击者类型和行动的组合，从收益矩阵 U 中获取相应的效用函数值，
                % 乘以攻击者类型的概率 P((t_A, t_D) 和攻击者类型的概率 F_star_A(t_A)，然后将结果加到 efficiency 变量上。
                efficiency(s_D_j) = efficiency(s_D_j) + U(s_A_i, s_D_j) * F_star(s_A_i);
            end       
        end
        efficiency_res = sum(efficiency .* F_strategy);
    elseif strcmp(player,'attacker')
        efficiency = zeros(1,nA);
        for s_A_i = 1:nA
            for s_D_j = 1:nD
                % 计算并累加效能期望值。对于每种攻击者类型和行动的组合，从收益矩阵 U 中获取相应的效用函数值，
                % 乘以攻击者类型的概率 P((t_A, t_D) 和攻击者类型的概率 F_star_A(t_A)，然后将结果加到 efficiency 变量上。
                efficiency(s_A_i) = efficiency(s_A_i) + U(s_A_i, s_D_j) * F_star(s_D_j);
            end
        end
        efficiency_res = sum(efficiency .* F_strategy);
    end
end

% 攻击者的目标函数
function J = objective_A(x)
    F_star_temp = x;
    
    % 计算不同策略的效能期望
    efficiency_A1_temp = calculate_efficiency_paper_Jiang(F_star_D, U, 1, 'attacker');
    efficiency_A2_temp = calculate_efficiency_paper_Jiang(F_star_D, U, 2, 'attacker');
    
    % 令攻击者的收益最小，取正号
    J = (efficiency_A1_temp*F_star_temp(1) + efficiency_A2_temp*F_star_temp(2));
end

% 防御者的目标函数
function J = objective_D(x)
    F_star_temp = x;
    
    % 计算不同策略的效能期望
    efficiency_D1_temp = calculate_efficiency_paper_Jiang(P_A, F_star_A, U, 1, 'attacker');
    efficiency_D2_temp = calculate_efficiency_paper_Jiang(P_A, F_star_A, U, 2, 'attacker');
    
    % 令防御者的收益最大，取负号
    J = -(efficiency_D1_temp*F_star_temp(1) + efficiency_D2_temp*F_star_temp(2));
end
end