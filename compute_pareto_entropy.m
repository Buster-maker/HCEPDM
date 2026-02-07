function entropy = compute_pareto_entropy(population)
    % population: 种群矩阵，每行代表一个解
    % 计算种群中每个维度的Pareto熵
    [num_individuals, num_dimensions] = size(population);
    entropy = zeros(1, num_dimensions);
    
    for d = 1:num_dimensions
        % 获取当前维度的所有解
        dimension_values = population(:, d);
        
        % 计算该维度的熵
        [counts, edges] = histcounts(dimension_values, 10); % 10个区间
        prob = counts / num_individuals;
        prob(prob == 0) = [];  % 移除零概率的项
        
        % 计算熵
        entropy(d) = -sum(prob .* log(prob)); 
    end
end
