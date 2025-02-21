clc;clear;close all;

% ���ġ�Near-optimal evalutation of network survivability under multi-stage attacks����֤
% �����������������������е����ݹ���
% �б�ʾ�����߲��ԣ��б�ʾ�����߲���
DOD_matrix = [
    0.749355, 2.979505, 2.99328, 2.995515, 2.996635;
    0.854925, 1.870855, 2.305035, 2.448675, 2.449975;
    0.99734, 1.632309, 1.872685, 2.0907, 2.148695;
    1.197075, 1.624254, 1.821465, 1.872965, 1.84314;
    1.496635, 1.798021, 1.997645, 2.13829, 1.49871
];

% ��ȡ�����߲��Ե������������������
num_defender_strategies = size(DOD_matrix, 1);
% ��ȡ�����߲��Ե������������������
num_attacker_strategies = size(DOD_matrix, 2);

% ���ڱ�Ƿ����ߵ�ÿ�������Ƿ�Ϊ���Ʋ��ԣ���ʼ����Ϊ false
dominant_defender_strategy = false(num_defender_strategies, 1);
% ���������ߵ�ÿһ������
for i = 1:num_defender_strategies
    % �ȼ��赱ǰ���� i �����Ʋ���
    is_dominant = true;
    % ����ǰ���� i ���������з����߲��Խ��бȽ�
    for j = 1:num_defender_strategies
        % ����������ıȽ�
        if i ~= j
            % ���ڹ����ߵ�ÿһ������
            for k = 1:num_attacker_strategies
                % �����ĳ�������߲��� k �£���ǰ���� i �� DOD ֵ���ڲ��� j �� DOD ֵ
                % ˵����ǰ���� i �������Ʋ���
                if DOD_matrix(i, k) > DOD_matrix(j, k)
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
    dominant_defender_strategy(i) = is_dominant;
end

% ���ڱ�ǹ����ߵ�ÿ�������Ƿ�Ϊ���Ʋ��ԣ���ʼ����Ϊ false
dominant_attacker_strategy = false(num_attacker_strategies, 1);
% ���������ߵ�ÿһ������
for i = 1:num_attacker_strategies
    % �ȼ��赱ǰ���� i �����Ʋ���
    is_dominant = true;
    % ����ǰ���� i ���������й����߲��Խ��бȽ�
    for j = 1:num_attacker_strategies
        % ����������ıȽ�
        if i ~= j
            % ���ڷ����ߵ�ÿһ������
            for k = 1:num_defender_strategies
                % �����ĳ�������߲��� k �£���ǰ���� i �� DOD ֵС�ڲ��� j �� DOD ֵ
                % ˵����ǰ���� i �������Ʋ���
                if DOD_matrix(k, i) < DOD_matrix(k, j)
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
    dominant_attacker_strategy(i) = is_dominant;
end

% ���û�����Ʋ��ԣ�ʹ��ԭ���������С���������
if ~any(dominant_defender_strategy) && ~any(dominant_attacker_strategy)
    % ��ʼ�������ߺͷ����ߵ����Ų�������
    attacker_opt_strategy = 1;
    defender_opt_strategy = 1;

    % ��ʼ����С���������Сֵ
    max_min_attacker = -Inf;
    min_max_defender = Inf;

    % Ѱ�ҹ����ߵ����Ų��ԣ������С���棩
    for i = 1:num_attacker_strategies
        min_value = min(DOD_matrix(:, i));
        if min_value > max_min_attacker
            max_min_attacker = min_value;
            attacker_opt_strategy = i;
        end
    end

    % Ѱ�ҷ����ߵ����Ų��ԣ���С�������ʧ��
    for j = 1:num_defender_strategies
        max_value = max(DOD_matrix(j, :));
        if max_value < min_max_defender
            min_max_defender = max_value;
            defender_opt_strategy = j;
        end
    end
end

% ��������ȷ���������Դ�������
attacker_strategies = {[0, 20], [5, 15], [10, 10], [15, 5], [20, 0]};
defender_strategies = {[0, 20], [5, 15], [10, 10], [15, 5], [20, 0]};
optimal_attacker_strategy = attacker_strategies{attacker_opt_strategy};
optimal_defender_strategy = defender_strategies{defender_opt_strategy};

disp(['�����ߵ����Ų���: ', num2str(optimal_attacker_strategy)]);
disp(['�����ߵ����Ų���: ', num2str(optimal_defender_strategy)]);