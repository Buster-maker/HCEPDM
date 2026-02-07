function chromosome = chaotic_initialize_variables(pop, M, V, min_range, max_range, item, tao, cntao, probID)
    % % Logistic映射参数
    % r = 3.9; % 混沌参数
    % x0 = rand; % 随机初始值
    % 
    % % 生成混沌序列
    % chaotic_sequence = zeros(pop * V, 1);
    % chaotic_sequence(1) = x0;
    % for i = 2:pop * V
    %     chaotic_sequence(i) = r * chaotic_sequence(i - 1) * (1 - chaotic_sequence(i - 1));
    % end
    % 
    % % 映射到变量范围
    % chromosome = zeros(pop, V + M);
    % for i = 1:pop
    %     for j = 1:V
    %         index = (i - 1) * V + j;
    %         chromosome(i, j) = min_range(j) + chaotic_sequence(index) * (max_range(j) - min_range(j));
    %     end
    % end
    % 高斯混沌映射参数
    a = 1.4;
    b = 0.3;

    % 生成高斯混沌序列
    chaotic_sequence = zeros(pop * V, 1);
    chaotic_sequence(1) = rand; % 随机初始值
    for i = 2:pop * V
        chaotic_sequence(i) = exp(-a * chaotic_sequence(i-1)^2) + b;
    end

    % 映射到变量范围
    chromosome = zeros(pop, V + M);
    for i = 1:pop
        for j = 1:V
            index = (i - 1) * V + j;
            chromosome(i, j) = min_range(j) + chaotic_sequence(index) * (max_range(j) - min_range(j));
        end
    end

