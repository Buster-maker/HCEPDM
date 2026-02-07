% 计算 Pareto 熵值
function [paretoEntropy, deltaEntropy] = calculateParetoEntropy(ParetoS)
    % 计算种群的熵值和熵变化
    [N, ~] = size(ParetoS);
    % 将数据划分为10个区间
    binCounts = histcounts(ParetoS, 10);
    % 计算每个区间的概率
    p = binCounts / N;
    % 计算熵值，值越高，解集的分布越均匀，多样性越好
    paretoEntropy = -sum(p .* log(p + eps));
    % 计算熵变化，即解集分布的最大概率与最小概率之差，反映了解集分布的离散程度
    deltaEntropy = max(p) - min(p);
end