clc;clear;close all;

% ���ġ���̬��Ҷ˹����������������ѡȡ��������֤
% �����������
P_D = [1,0];        % �����߶ԡ������ߡ����͵���������
P_A = [0.5,0.5];    % �����߶ԡ������ߡ����͵���������

% ��ϲ���
F_star_D = [0.5, 0.5];
F_star_A = [0.2, 0.8];

% �����������ÿ��Ԫ�� (m_A_i, s_D_j, t_A) ��Ӧһ������ֵ
% ÿ������������У��� 1 ��Ϊ�����ߵĲ������棬�� 2 ����Ϊ�����ߵĲ�������
U = cat(3, ...
    {[40,50], [-10,0]; [0,300], [0,300]}, ...   % ����������1
    {[30,80], [-10,100]; [0,400], [0,400]});    % ����������2

% ��������Ч��
efficiency_D1 = calculate_efficiency(P_D, F_star_D, U, 1, 'defender');
efficiency_D2 = calculate_efficiency(P_D, F_star_D, U, 2, 'defender');

fprintf('E_s_D1: %.4f\n', efficiency_D1);
fprintf('E_s_D2: %.4f\n', efficiency_D2);

% ��������Ч��
efficiency_A1 = calculate_efficiency(P_A, F_star_A, U, 1, 'attacker');
efficiency_A2 = calculate_efficiency(P_A, F_star_A, U, 2, 'attacker');

fprintf('E_s_A1: %.4f\n', efficiency_A1);
fprintf('E_s_A2: %.4f\n', efficiency_A2);

function efficiency = calculate_efficiency(P, F_star, U, strategy, player)
% ˵����������������Ч��
% ���㹫ʽ������̬��Ҷ˹����������������ѡȡ������ʽ(1)
% ����˵����
% ���룺
% P ���� �ҷ��Եз����͵��������
% F_star ���� �з�����ѡȡ����
% U ���� �����������
% strategt ���� �ҷ�ѡȡ�ĵ�i������
% player ���� �ҷ���ݣ�������defender���͡�attacker��
% �����
% efficiency ���� ��������Ч��
%--------------------------------------------------------------------------
    efficiency = 0;
    if strcmp(player,'defender')
        nD = length(strategy);  % Defender�Ĳ�����
        nA = length(F_star);    % Attacker�Ĳ�����
        player_mode = 1;        % �� U ��ʹ�ã�ȡ�������ĵ� 1 λ��
        for t_A = 1:nA          % �������п��ܵĹ��������� t_A
            for m_A_i = 1:nD    % �������п��ܵĹ����߲��� m_A_i
                % ���㲢�ۼ�Ч������ֵ������ÿ�ֹ��������ͺ��ж�����ϣ���������� U �л�ȡ��Ӧ��Ч�ú���ֵ��
                % ���Թ��������͵ĸ��� P((t_A, t_D) �͹��������͵ĸ��� F_star_A(t_A)��Ȼ�󽫽���ӵ� efficiency �����ϡ�
                efficiency = efficiency + P(t_A) * U{m_A_i, strategy, t_A}(player_mode) * F_star(m_A_i);
            end
        end                
    elseif strcmp(player,'attacker')
        nD = length(strategy);  % Defender�Ĳ�����
        nA = length(F_star);    % Attacker�Ĳ�����        
        player_mode = 2;        % �� U ��ʹ�ã�ȡ�������ĵ� 2 λ��
        for t_A = 1:nA          % �������п��ܵĹ��������� t_A
            for s_D_j = 1:nD    % �������п��ܵķ����߲��� m_A_i
                % **** ����Ҫ�� **** ��ΪF_starҪ�÷����ߵģ�״̬��Done��ʱ�䣺2025/02/20��
                % ���㲢�ۼ�Ч������ֵ������ÿ�ֹ��������ͺ��ж�����ϣ���������� U �л�ȡ��Ӧ��Ч�ú���ֵ��
                % ���Թ��������͵ĸ��� P((t_A, t_D) �͹��������͵ĸ��� F_star_A(t_A)��Ȼ�󽫽���ӵ� efficiency �����ϡ�
                efficiency = efficiency + P(1) * U{strategy, s_D_j, t_A}(player_mode) * F_star(s_D_j);
            end
        end
    end
end
