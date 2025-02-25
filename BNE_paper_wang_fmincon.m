clc;clear;

%% 主函数
% 说明：复现论文《静态贝叶斯博弈主动防御策略选取方法》——fimcon()方法
% 定义概率，攻击者是冒险型的概率
p = 0.6;

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
U = cat(3,...
    high_benefit_matrix,...     % 攻击者为高技术情况
    low_benefit_matrix);        % 攻击者为低技术情况

% 获取 D 的策略数量（矩阵的行数）
nD = size(high_benefit_matrix, 1);
% 获取 A 的策略数量（矩阵的列数）
nA = size(high_benefit_matrix, 2);

% 计算期望博弈收益定义变量
P_A = [p,1-p];      % 防御者对【攻击者】类型的先验信念
P_D = [1,0];        % 攻击者对【防御者】类型的先验信念
num_type_A = 2;     % 攻击者类型数
num_type_D = 1;     % 防御者类型数
num_type = [num_type_A, num_type_D];    % 第1个数是攻击者类型数，第2个数是防御者类型数

% 初始化混合策略
F_star_D = ones(1,nD)/nD;
F_star_A_high = ones(1,nA)/nA;
F_star_A_low = ones(1,nA)/nA;
F_star_A = [F_star_A_high;F_star_A_low];

%% 优化求解纳什均衡
max_iter = 1000;
tol = 1e-2;

% 优化
for iter = 1:max_iter 
    F_star_D_temp = optimize_F_D(num_type, P_A, F_star_D, F_star_A, U, 'defender', 'none');
    F_star_A_high_temp = optimize_F_A_high(num_type, P_D, F_star_D, F_star_A_high, U, 'attacker', 'high');
    F_star_A_low_temp = optimize_F_A_low(num_type, P_D, F_star_D, F_star_A_low, U, 'attacker', 'low');
    
    if ( norm(F_star_A_high_temp-F_star_A(1,:)) )< tol && ( norm(F_star_A_low_temp-F_star_A(2,:)) )< tol && (norm(F_star_D_temp-F_star_D))< tol
            F_star_D = F_star_D_temp;
            F_star_A_high = F_star_A_high_temp;
            F_star_A_low = F_star_A_low_temp;
            F_star_A = [F_star_A_high_temp;F_star_A_low_temp];
        break;
    end
    % 更新变量
    F_star_D = F_star_D_temp;
    F_star_A_high = F_star_A_high_temp;
    F_star_A_low = F_star_A_low_temp;
    F_star_A = [F_star_A_high_temp;F_star_A_low_temp];
end

% 返回策略函数
disp(F_star_A_high_temp);
disp(F_star_A_low_temp);
disp(F_star_D);

%% 子函数部分
function sum_U = utility(num_type,P,F_D,F_A,U,player,mode)
    num_type_A = num_type(1);   % 攻击者类型数
    num_type_D = num_type(2);   % 防御者类型数
    switch player
        case 'defender'
            nD = length(F_D);
            sum_U_temp = zeros(num_type_A,1);
            player_mode = 1;        % 在 U 里使用，取收益矩阵的第 1 位数
            benefit_matrix_high = cell2mat(U(:,:,1));
            benefit_matrix_high = benefit_matrix_high(:,player_mode:2:nD*2);
            
            benefit_matrix_low = cell2mat(U(:,:,2));
            benefit_matrix_low = benefit_matrix_low(:,player_mode:2:nD*2);
            
            sum_U_temp(1) = P(1) * F_D * benefit_matrix_high * F_A(1,:)';
            sum_U_temp(2) = P(2) * F_D * benefit_matrix_low * F_A(2,:)';
            sum_U = sum(sum_U_temp);
        case 'attacker'
            player_mode = 2;        % 在 U 里使用，取收益矩阵的第 2 位数
            nA = length(F_A);
            switch mode
                case 'high'
                    benefit_matrix_high = cell2mat(U(:,:,1));
                    benefit_matrix_high = benefit_matrix_high(:,player_mode:2:nA*2);
                    sum_U = P(num_type_D) * ( F_D * benefit_matrix_high * F_A' );
                case 'low'
                    benefit_matrix_low = cell2mat(U(:,:,2));
                    benefit_matrix_low = benefit_matrix_low(:,player_mode:2:nA*2);
                    sum_U = P(num_type_D) * ( F_D * benefit_matrix_low * F_A' );
            end
            
    end
            
end

% 定义优化函数
function F_opt_A_low = optimize_F_A_low(num_type,P,F_D,F_A,U,player,mode)
    options = optimoptions('fmincon', 'Display', 'off');
    nAA = length(F_A);
    F_opt_A_low = fmincon(@objective_F_A, F_A, [], [], ones(1,nAA), 1, zeros(1,nAA), ones(1,nAA), [], options);
    function obj = objective_F_A(F_A)
        obj = -utility(num_type,P,F_D,F_A,U,player,mode);
    end
end


function F_opt_A_high = optimize_F_A_high(num_type,P,F_D,F_A,U,player,mode)
    options = optimoptions('fmincon', 'Display', 'off');
    nAA = length(F_A);
    F_opt_A_high = fmincon(@objective_F_A, F_A, [], [], ones(1,nAA), 1, zeros(1,nAA), ones(1,nAA), [], options);
    function obj = objective_F_A(F_A)
        obj = -utility(num_type,P,F_D,F_A,U,player,mode);
    end
end

function F_opt_D = optimize_F_D(num_type,P,F_D,F_A,U,player,mode)
    options = optimoptions('fmincon', 'Display', 'off');
    nDD = length(F_D);
    F_opt_D = fmincon(@objective_F_D, F_D, [], [], ones(1,nDD), 1, zeros(1,nDD), ones(1,nDD), [], options);
    function obj = objective_F_D(F_D)
        obj = -utility(num_type,P,F_D,F_A,U,player,mode);
    end
end
