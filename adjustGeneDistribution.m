% 调整基因值的均匀分布函数
function chromosome = adjustGeneDistribution(chromosome, probID, min_range, max_range)
    CAJ = CaJ_erwei(probID); % 获取需要调整的维度和其他参数
    [pop, V] = size(chromosome); % 获取种群大小和决策变量数

    % 获取最大值和最小值
    [~, index1] = max(chromosome);
    [~, index2] = min(chromosome);
    max1 = chromosome(index1(CAJ(1, 2)), CAJ(1, 2));
    min1 = chromosome(index2(CAJ(1, 2)), CAJ(1, 2));
    
    % 对指定维度的基因值进行均匀分布调整
    for w = 1:pop
        BY = min1 + ((max1 - min1) / pop) * w;
        chromosome(w, CAJ(1, 2)) = BY;
    end
    
    % 确保基因值在范围内
    chromosome(:, CAJ(1, 2)) = max(min(chromosome(:, CAJ(1, 2)), max_range(CAJ(1, 2))), min_range(CAJ(1, 2)));
end