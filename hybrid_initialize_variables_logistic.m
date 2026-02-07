function chromosome = hybrid_initialize_variables_logistic(pop, M, V, min_range, max_range, item, tao, cntao, probID)
    % 初始化参数
    ratioChaos = 0.6; % Logistic混沌比例
    ratioRandom = 0.4; % 均匀分布比例
    
    popChaos = round(pop * ratioChaos);
    popRandom = pop - popChaos;

    % Logistic混沌映射初始化
    r = 3.9;
    x = rand;
    chaotic_sequence = zeros(popChaos * V, 1);
    chaotic_sequence(1) = x;
    for i = 2:popChaos * V
        chaotic_sequence(i) = r * chaotic_sequence(i-1) * (1 - chaotic_sequence(i-1));
    end
    chaos_samples = zeros(popChaos, V);
    for i = 1:popChaos
        for j = 1:V
            index = (i - 1) * V + j;
            chaos_samples(i, j) = min_range(j) + chaotic_sequence(index) * (max_range(j) - min_range(j));
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
    decision_variables = [chaos_samples; random_samples];
    
    % 计算目标函数值
    chromosome = zeros(pop, V + M);
    for i = 1:pop
        chromosome(i, 1:V) = decision_variables(i, :);
        chromosome(i, V + 1:V + M) = evaluate_objective(chromosome(i, 1:V), M, V, item, tao, cntao, probID);
    end
end

function objective_values = evaluate_objective(decision_variables, M, V, item, tao, cntao, probID)
    % 这里是计算目标函数值的代码，根据您的具体问题和目标函数来实现
    % 例如：
    objective_values = zeros(1, M);
    for i = 1:M
        % 假设目标函数是所有决策变量的平方和
        objective_values(i) = sum(decision_variables.^2);
    end
end
