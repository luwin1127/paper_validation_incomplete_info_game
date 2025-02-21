function Bayesian_equilibrium_selection_paper_Jiang
% ��֤���ġ����ڹ�������ģ�͵����簲ȫ����������������������ΰ
clc;clear;close all;

%% ��ʼ��������
% ���ġ���̬��Ҷ˹����������������ѡȡ��������֤
% �����������
P_D = [1,0];        % �����߶ԡ������ߡ����͵���������
P_A = [0.9,0.1];    % �����߶ԡ������ߡ����͵���������
type_P_D = 1;       % ������������
type_P_A = 2;       % ������������

% ��������
nD = 6; % �����߲�������
nA = 3; % �����߲�������

% ��ϲ���
% F_star_D = ones(1,nD)*(1/nD);
% F_star_A = ones(1,nA)*(1/nA);
% ����Ľ�������Ž��
F_star_D = [28/53, 0, 0, 0, 0, 25/36];
F_star_A = [67/159, 0, 92/159];

% �����������ÿ��Ԫ�� (s_A, s_D) ��Ӧһ������ֵ
U = [810 1510 1670 2470 2650 2300
     810 2470 1670 2470 2470 2300
    2310 2310 2470 1670  970 2300];

% ��������Ч��
efficiency_D = calculate_efficiency_paper_Jiang(F_star_A, U, F_star_D, 'defender');
fprintf('E_s_D: %.4f\n', efficiency_D);

% ��������Ч��
efficiency_A = calculate_efficiency_paper_Jiang(F_star_D, U, F_star_A, 'attacker');
fprintf('E_s_A1: %.4f\n', efficiency_A);

%% �㷨��ʼ
% ѭ������
% tol = 1;
options = optimset('Display','off',...  % off/iter/notify-detailed
                   'Algorithm','sqp',...% �㷨ѡ��
                   'TolCon',1e-6,...    % Լ��Υ���ȶ��ݲ����������Ĭ��ֵΪ 1e-6
                   'TolX',1e-6,...      % ���������� x ����ֹ�ݲSQP�㷨��Ĭ��ֵΪ 1e-6
                   'TolFun',1e-6);      % һ�������Ե���ֹ�ݲ����������Ĭ��ֵΪ 1e-6             
lb = [0;0];
ub = [1;1];
Aeq = [1 1];
beq = 1;

% ����ѡȡ�㷨
x_D = fmincon(@objective_D,F_star_D,[],[],Aeq,beq,lb,ub,[],options);
F_star_D = x_D;
disp(F_star_D);

x_A = fmincon(@objective_A,F_star_A,[],[],Aeq,beq,lb,ub,[],options);
F_star_A = x_A;
disp(F_star_A);

%% �����ʾ
% ��������Ч��
efficiency_D1 = calculate_efficiency_paper_Jiang(F_star_A, U, 1, 'defender');
efficiency_D2 = calculate_efficiency_paper_Jiang(F_star_A, U, 2, 'defender');

fprintf('E_s_D1: %.4f\n', efficiency_D1);
fprintf('E_s_D2: %.4f\n', efficiency_D2);

% ��������Ч��
efficiency_A1 = calculate_efficiency_paper_Jiang(F_star_D, U, 1, 'attacker');
efficiency_A2 = calculate_efficiency_paper_Jiang(F_star_D, U, 2, 'attacker');

fprintf('E_s_A1: %.4f\n', efficiency_A1);
fprintf('E_s_A2: %.4f\n', efficiency_A2);
%% �Ӻ�������
% ����Ч������
function efficiency_res = calculate_efficiency_paper_Jiang(F_star, U, F_strategy, player)
    if strcmp(player,'defender')
        efficiency = zeros(1,nD);
        for s_D_j = 1:nD
            for s_A_i = 1:nA % �������п��ܵĹ����߲��� m_A_i
                % ���㲢�ۼ�Ч������ֵ������ÿ�ֹ��������ͺ��ж�����ϣ���������� U �л�ȡ��Ӧ��Ч�ú���ֵ��
                % ���Թ��������͵ĸ��� P((t_A, t_D) �͹��������͵ĸ��� F_star_A(t_A)��Ȼ�󽫽���ӵ� efficiency �����ϡ�
                efficiency(s_D_j) = efficiency(s_D_j) + U(s_A_i, s_D_j) * F_star(s_A_i);
            end       
        end
        efficiency_res = sum(efficiency .* F_strategy);
    elseif strcmp(player,'attacker')
        efficiency = zeros(1,nA);
        for s_A_i = 1:nA
            for s_D_j = 1:nD
                % ���㲢�ۼ�Ч������ֵ������ÿ�ֹ��������ͺ��ж�����ϣ���������� U �л�ȡ��Ӧ��Ч�ú���ֵ��
                % ���Թ��������͵ĸ��� P((t_A, t_D) �͹��������͵ĸ��� F_star_A(t_A)��Ȼ�󽫽���ӵ� efficiency �����ϡ�
                efficiency(s_A_i) = efficiency(s_A_i) + U(s_A_i, s_D_j) * F_star(s_D_j);
            end
        end
        efficiency_res = sum(efficiency .* F_strategy);
    end
end

% �����ߵ�Ŀ�꺯��
function J = objective_A(x)
    F_star_temp = x;
    
    % ���㲻ͬ���Ե�Ч������
    efficiency_A1_temp = calculate_efficiency_paper_Jiang(F_star_D, U, 1, 'attacker');
    efficiency_A2_temp = calculate_efficiency_paper_Jiang(F_star_D, U, 2, 'attacker');
    
    % ����ߵ�������С��ȡ����
    J = (efficiency_A1_temp*F_star_temp(1) + efficiency_A2_temp*F_star_temp(2));
end

% �����ߵ�Ŀ�꺯��
function J = objective_D(x)
    F_star_temp = x;
    
    % ���㲻ͬ���Ե�Ч������
    efficiency_D1_temp = calculate_efficiency_paper_Jiang(P_A, F_star_A, U, 1, 'attacker');
    efficiency_D2_temp = calculate_efficiency_paper_Jiang(P_A, F_star_A, U, 2, 'attacker');
    
    % ������ߵ��������ȡ����
    J = -(efficiency_D1_temp*F_star_temp(1) + efficiency_D2_temp*F_star_temp(2));
end
end