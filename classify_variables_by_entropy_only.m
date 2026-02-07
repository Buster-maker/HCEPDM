function [classification, entropy_std] = classify_variables_by_entropy_only(PF_history)
% 基于帕累托解的熵波动对决策变量分类（仅基于分布，不依赖目标函数）

[T, N, V] = size(PF_history);
entropy_values = zeros(T, V);

for t = 1:T
    data = reshape(PF_history(t, 1:N, :), [N, V]);
    for v = 1:V
        col = data(:, v);
        if all(abs(col - col(1)) < 1e-8)
            entropy_values(t, v) = 0;
        else
            [~, edges] = histcounts(col, 10);
            counts = histcounts(col, edges, 'Normalization', 'probability');
            counts(counts == 0) = [];
            entropy_values(t, v) = -sum(counts .* log(counts + 1e-10));
        end
    end
end

entropy_std = std(entropy_values, 0, 1);
mean_entropy = mean(entropy_std);

classification = strings(1, V);
for v = 1:V
    if entropy_std(v) >= mean_entropy
        classification(v) = "多样性";
    else
        classification(v) = "收敛性";
    end
end
end
