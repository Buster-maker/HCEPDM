function [classification, entropy_std, sensitivity] = classify_variables_by_entropy_sensitivity(PF_history, obj_fun, min_range, max_range, num_samples, epsilon)
% 基于帕累托前沿熵响应和敏感性分析的决策变量分类方法
% 输入:
%   PF_history: T × N × V 的3D矩阵，记录T代的帕累托解
%   obj_fun: 匿名函数 @(x) evaluate_objective(...)
%   min_range, max_range: 决策变量范围
%   num_samples: 每变量敏感性扰动次数
%   epsilon: 扰动步长
% 输出:
%   classification: 每个变量的类别（多样性变量/收敛性变量/中性变量）
%   entropy_std: 每个变量在多代中的熵波动（标准差）
%   sensitivity: 每个变量对目标函数的扰动响应均值

[T, N, V] = size(PF_history);

entropy_values = zeros(T, V);

for t = 1:T
    data = reshape(PF_history(t, 1:N, :), [N, V]);  % N × V 保证不降维
    for v = 1:V
        col = data(:, v);

        % 若变量为常值，熵为0
        if all(abs(col - col(1)) < 1e-8)
            entropy_values(t, v) = 0;
            continue;
        end

        try
            [~, edges] = histcounts(col, 10);
            counts = histcounts(col, edges, 'Normalization', 'probability');
            counts(counts == 0) = [];
            entropy_values(t, v) = -sum(counts .* log(counts + 1e-10));
        catch
            entropy_values(t, v) = 0;
        end
    end
end

% 计算标准差（熵波动）
entropy_std = std(entropy_values, 0, 1);  % 1 × V

%% --- 敏感性分析 ---
sensitivity = zeros(1, V);
base_sample = (min_range + max_range) / 2;

for v = 1:V
    diffs = zeros(1, num_samples);
    for i = 1:num_samples
        x1 = base_sample;
        x2 = base_sample;
        x1(v) = min(max(x1(v) - epsilon, min_range(v)), max_range(v));
        x2(v) = min(max(x2(v) + epsilon, min_range(v)), max_range(v));
        try
            y1 = obj_fun(x1);
            y2 = obj_fun(x2);
            diffs(i) = norm(y2 - y1);
        catch
            diffs(i) = 0;
        end
    end
    sensitivity(v) = mean(diffs);
end

%% --- 分类规则 ---
classification = strings(1, V);
entropy_z = (entropy_std - mean(entropy_std)) ./ std(entropy_std);
sens_z    = (sensitivity - mean(sensitivity)) ./ std(sensitivity);

for v = 1:V
    if entropy_z(v) >= -0.25 && sens_z(v) <= 0.25
        classification(v) = "多样性变量";
    else
        classification(v) = "收敛性变量";
    end
end


end
