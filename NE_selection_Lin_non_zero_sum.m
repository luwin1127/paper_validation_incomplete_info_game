clc;clear;close all;

% ���ġ�Near-optimal evalutation of network survivability under multi-stage attacks����չΪ����ȫ��Ϣ��̬����
% ������ʣ���������ǿ���͵ĸ���
p = 0.6;

% ����ǿ���͹������µ��������ʹ��Ԫ������洢
DOD_matrix_strong = {
    [2, 8], [1, 9], [0.5, 9.5], [0.2, 9.8], [0.1, 9.9];
    [2.5, 7.5], [2, 8], [1.5, 8.5], [1.2, 8.8], [1.1, 8.9];
    [3, 7], [2.5, 7.5], [2, 8], [1.8, 8.2], [1.7, 8.3];
    [3.5, 6.5], [3, 7], [2.5, 7.5], [2.2, 7.8], [2.1, 7.9];
    [4, 6], [3.5, 6.5], [3, 7], [2.8, 7.2], [2.7, 7.3]
};

% ���������͹������µ��������ʹ��Ԫ������洢
DOD_matrix_weak = {
    [3, 7], [2, 8], [1.5, 8.5], [1.2, 8.8], [1.1, 8.9];
    [3.5, 6.5], [3, 7], [2.5, 7.5], [2.2, 7.8], [2.1, 7.9];
    [4, 6], [3.5, 6.5], [3, 7], [2.8, 7.2], [2.7, 7.3];
    [4.5, 5.5], [4, 6], [3.5, 6.5], [3.2, 6.8], [3.1, 6.9];
    [5, 5], [4.5, 5.5], [4, 6], [3.8, 6.2], [3.7, 6.3]
};

% ��ȡ�����߲��Ե������������������
num_defender_strategies = size(DOD_matrix_strong, 1);
% ��ȡ�����߲��Ե������������������
num_attacker_strategies = size(DOD_matrix_strong, 2);

% ��ʼ������������󣬷ֱ�洢�����ߺ͹����ߵ���������
expected_defender_matrix = zeros(num_defender_strategies, num_attacker_strategies);
expected_attacker_matrix = zeros(num_defender_strategies, num_attacker_strategies);

% ���������������
for i = 1:num_defender_strategies
    for j = 1:num_attacker_strategies
        expected_defender_matrix(i, j) = p * DOD_matrix_strong{i, j}(1) + (1 - p) * DOD_matrix_weak{i, j}(1);
        expected_attacker_matrix(i, j) = p * DOD_matrix_strong{i, j}(2) + (1 - p) * DOD_matrix_weak{i, j}(2);
    end
end

% ������Ʋ��Բ��� - ������
% ��ʼ��һ���������������ڱ�Ƿ����ߵ�ÿ�������Ƿ�Ϊ���Ʋ��ԣ���ʼ����Ϊ false
dominant_defender_strategy = false(num_defender_strategies, 1);
% ���ѭ�������������ߵ�ÿһ������
for i = 1:num_defender_strategies
    % �ȼ��赱ǰ���� i �����Ʋ���
    is_dominant = true;
    % �м��ѭ��������ǰ���� i ���������з����߲��Խ��бȽ�
    for j = 1:num_defender_strategies
        % ����������ıȽ�
        if i ~= j
            % �ڲ�ѭ�������ڹ����ߵ�ÿһ������
            for k = 1:num_attacker_strategies
                % �����ĳ�������߲��� k �£���ǰ���� i ����������С�ڲ��� j ����������
                % ˵����ǰ���� i �������Ʋ���
                if expected_defender_matrix(i, k) < expected_defender_matrix(j, k)
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

% ������Ʋ��Բ��� - ������
% ��ʼ��һ���������������ڱ�ǹ����ߵ�ÿ�������Ƿ�Ϊ���Ʋ��ԣ���ʼ����Ϊ false
dominant_attacker_strategy = false(num_attacker_strategies, 1);
% ���ѭ�������������ߵ�ÿһ������
for i = 1:num_attacker_strategies
    % �ȼ��赱ǰ���� i �����Ʋ���
    is_dominant = true;
    % �м��ѭ��������ǰ���� i ���������й����߲��Խ��бȽ�
    for j = 1:num_attacker_strategies
        % ����������ıȽ�
        if i ~= j
            % �ڲ�ѭ�������ڷ����ߵ�ÿһ������
            for k = 1:num_defender_strategies
                % �����ĳ�������߲��� k �£���ǰ���� i ����������С�ڲ��� j ����������
                % ˵����ǰ���� i �������Ʋ���
                if expected_attacker_matrix(k, i) < expected_attacker_matrix(k, j)
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

% ���û�����Ʋ��ԣ�ʹ�������������������
if ~any(dominant_defender_strategy) && ~any(dominant_attacker_strategy)
    % ��ʼ�������ߺͷ����ߵ����Ų�������
    attacker_opt_strategy = 1;
    defender_opt_strategy = 1;

    % ��ʼ�������������
    max_attacker_income = -Inf;
    max_defender_income = -Inf;

    % Ѱ�ҹ����ߵ����Ų��ԣ�����������棩
    for i = 1:num_attacker_strategies
        total_income = sum(expected_attacker_matrix(:, i));
        if total_income > max_attacker_income
            max_attacker_income = total_income;
            attacker_opt_strategy = i;
        end
    end

    % Ѱ�ҷ����ߵ����Ų��ԣ�����������棩
    for j = 1:num_defender_strategies
        total_income = sum(expected_defender_matrix(j, :));
        if total_income > max_defender_income
            max_defender_income = total_income;
            defender_opt_strategy = j;
        end
    end
else
    % ��������Ʋ��ԣ�ֱ�Ӹ������Ʋ���ȷ�����Ų���
    if any(dominant_defender_strategy)
        defender_opt_strategy = find(dominant_defender_strategy, 1);
    end
    if any(dominant_attacker_strategy)
        attacker_opt_strategy = find(dominant_attacker_strategy, 1);
    end
end

% ��������ȷ���������Դ�������
attacker_strategies = {[0, 20], [5, 15], [10, 10], [15, 5], [20, 0]};
defender_strategies = {[0, 20], [5, 15], [10, 10], [15, 5], [20, 0]};
optimal_attacker_strategy = attacker_strategies{attacker_opt_strategy};
optimal_defender_strategy = defender_strategies{defender_opt_strategy};

disp(['�����ߵ����Ų���: ', num2str(optimal_attacker_strategy)]);
disp(['�����ߵ����Ų���: ', num2str(optimal_defender_strategy)]);