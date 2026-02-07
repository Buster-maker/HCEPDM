function chromosome = disturb_population_by_variable_type(chromosome, classification, min_range, max_range, center, rate)
% 对不同变量维度采用不同扰动：多样性变量 → 均匀化扰动，收敛性变量 → 局部扰动

[pop, dim] = size(chromosome);
M = dim - length(classification); % 决策变量数V
V = length(classification);

for i = 1:pop
    for v = 1:V
        if rand < rate
            if classification(v) == "多样性"
                % 在全局范围均匀扰动
                chromosome(i, M+v) = min_range(v) + rand * (max_range(v) - min_range(v));
            elseif classification(v) == "收敛性"
                % 在中心附近进行局部扰动
                local_eps = 0.05 * (max_range(v) - min_range(v));
                chromosome(i, M+v) = center(v) + local_eps * (2*rand - 1);
                chromosome(i, M+v) = max(min(chromosome(i, M+v), max_range(v)), min_range(v));
            end
        end
    end
end
end
