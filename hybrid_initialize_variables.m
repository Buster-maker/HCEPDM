function chromosome = hybrid_initialize_variables(pop, M, V, min_range, max_range, item, tao, cntao, probID)
    % 初始化参数
    ratioLHS = 0.4; % LHS比例
    ratioChaos = 0.4; % 高斯混沌比例
    ratioRandom = 0.2; % 均匀分布比例
    
    popLHS = round(pop * ratioLHS);
    popChaos = round(pop * ratioChaos);
    popRandom = pop - popLHS - popChaos;

    % LHS初始化
    LHS_samples = lhsdesign(popLHS, V) .* (max_range - min_range) + min_range;

    % 高斯混沌映射初始化
    a = 1.4;
    b = 0.3;
    x = rand;
    y = rand;
    chaotic_sequence_x = zeros(popChaos * V, 1);
    chaotic_sequence_y = zeros(popChaos * V, 1);
    chaotic_sequence_x(1) = x;
    chaotic_sequence_y(1) = y;
    for i = 2:popChaos * V
        chaotic_sequence_x(i) = 1 - a * chaotic_sequence_x(i-1)^2 + chaotic_sequence_y(i-1);
        chaotic_sequence_y(i) = b * chaotic_sequence_x(i-1);
    end
    chaos_samples = zeros(popChaos, V);
    for i = 1:popChaos
        for j = 1:V
            index = (i - 1) * V + j;
            chaos_samples(i, j) = min_range(j) + chaotic_sequence_x(index) * (max_range(j) - min_range(j));
            chaos_samples(i, j) = max(min(chaos_samples(i, j), max_range(j)), min_range(j));
        end
    end

    % 均匀分布初始化
    random_samples = zeros(popRandom, V);
    for i = 1:popRandom
        for j = 1:V
            random_samples(i, j) = min_range(j) + (max_range(j) - min_range(j)) * rand();
        end
    end

    % 合并所有样本
    decision_variables = [LHS_samples; chaos_samples; random_samples];
    
    % 计算目标函数值
    chromosome = zeros(pop, V + M);
    for i = 1:pop
        chromosome(i, 1:V) = decision_variables(i, :);
        chromosome(i, V + 1:V + M) = evaluate_objective(chromosome(i, 1:V), M, V, item, tao, cntao, probID);
    end
end
