function chromosome = initialize_variables_with_chaos(pop, M, V, min_range, max_range, item, tao, cntao, probID)
    % 拉丁超立方体采样生成初始种群
    LHS_samples = lhsdesign(pop, V) .* (max_range - min_range) + min_range;
    
    % 高斯混沌映射参数
    a = 1.2;
    b = 0.3;
    x = rand; % 初始值
    y = rand; % 初始值

    % 生成混沌序列
    chaotic_sequence_x = zeros(pop * V, 1);
    chaotic_sequence_y = zeros(pop * V, 1);
    chaotic_sequence_x(1) = x;
    chaotic_sequence_y(1) = y;
    for i = 2:pop * V
        chaotic_sequence_x(i) = 1 - a * chaotic_sequence_x(i-1)^2 + chaotic_sequence_y(i-1);
        chaotic_sequence_y(i) = b * chaotic_sequence_x(i-1);
    end

    % 映射到变量范围
    chromosome = zeros(pop, V + M);
    for i = 1:pop
        for j = 1:V
            index = (i - 1) * V + j;
            % 将拉丁超立方体采样和混沌映射结合起来
            chromosome(i, j) = LHS_samples(i, j) + chaotic_sequence_x(index) * (max_range(j) - min_range(j));
            % 保证染色体在决策空间内
            chromosome(i, j) = max(min(chromosome(i, j), max_range(j)), min_range(j));
        end
    end
end

