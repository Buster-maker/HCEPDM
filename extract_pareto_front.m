function pareto_set = extract_pareto_front(chromosome, M, V)
    % 提取种群中的Pareto前沿解
    pareto_front = chromosome(:, 1:M); % 假设前M列为目标值
    pareto_set = chromosome(pareto_front(:, end) == 1, 1:V); % 假设最后一列为支配标志
end