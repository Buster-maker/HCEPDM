% 根据熵值修改种群
function chromosome = modifyPopulationBasedOnEntropy(chromosome, paretoEntropy, min_range, max_range, threshold_low)
    [pop, V] = size(chromosome);
    if paretoEntropy < threshold_low
        % 当熵值较低时，随机初始化部分个体以增加多样性
        num_random = round(pop * 0.1); % 例如随机初始化10%的个体
        random_indices = randperm(pop, num_random);
        for i = 1:num_random
            idx = random_indices(i);
            for j = 1:V
                chromosome(idx, j) = min_range(j) + (max_range(j) - min_range(j)) * rand();
            end
        end
    end
end
