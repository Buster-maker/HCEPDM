function f = initialize_variables_with_lhs_and_chaos_feedback(N, M, V, min_range, max_range, item, tao, cntao, probID, prev_entropy)
% 自适应熵反馈驱动的混合初始化方法
% 基于上一代Pareto熵动态分配子群混沌映射方式

% ---------- 参数 ----------
lhs_sample = lhsdesign(N, V); % 拉丁超立方初始分布
initial_population = zeros(N, V);
for j = 1:V
    initial_population(:, j) = min_range(j) + (max_range(j) - min_range(j)) .* lhs_sample(:, j);
end

% ---------- 子群划分 ----------
num_subgroups = 3;
subgroup_size = floor(N / num_subgroups);
f = zeros(N, M + V);

% ---------- 动态扰动幅度 ----------
base_perturbation = 0.05;
max_perturbation = 0.2;

% ---------- 根据熵值动态调整策略 ----------
% 设定熵阈值
low_entropy = 0.5;
high_entropy = 1.5;

% 熵低：使用高强度扰动（Tent优先）
% 熵高：使用稳定扰动（Gaussian优先）
if prev_entropy <= low_entropy
    chaos_order = {@tent_chaos, @logistic_chaos, @gaussian_chaos}; % 探索优先
elseif prev_entropy >= high_entropy
    chaos_order = {@gaussian_chaos, @logistic_chaos, @tent_chaos}; % 稳定收敛
else
    chaos_order = {@logistic_chaos, @gaussian_chaos, @tent_chaos}; % 均衡
end

% ---------- 执行扰动 ----------
for g = 1:num_subgroups
    start_idx = (g - 1) * subgroup_size + 1;
    end_idx = (g == num_subgroups) * N + (g ~= num_subgroups) * (g * subgroup_size);
    subgroup = initial_population(start_idx:end_idx, :);

    chaotic_seq = zeros(size(subgroup));
    for i = 1:size(subgroup, 1)
        chaotic_seq(i, :) = chaos_order{g}(subgroup(i, :));
    end

    % 动态扰动比例
    dynamic_perturbation = base_perturbation + ...
        (max_perturbation - base_perturbation) * rand(size(subgroup, 1), 1);
    dynamic_perturbation = repmat(dynamic_perturbation, 1, V);

    % 添加扰动并加随机微噪声
    perturbed = subgroup + dynamic_perturbation .* chaotic_seq + 0.01 * randn(size(subgroup));
    perturbed = max(min(perturbed, max_range), min_range); % 边界处理
    f(start_idx:end_idx, 1:V) = perturbed;
end

% ---------- 目标函数计算 ----------
for i = 1:N
    f(i, V + 1:end) = evaluate_objective(f(i, 1:V), M, V, item, tao, cntao, probID);
end
end

% Logistic混沌映射
function chaotic = logistic_chaos(values, r)
if nargin < 2
    r = 3.9;
end
chaotic = r .* values .* (1 - values);
end

% 高斯混沌映射
function chaotic = gaussian_chaos(values, a)
if nargin < 2
    a = 6;
end
chaotic = exp(-a .* values.^2) + 0.5;
end

% Tent混沌映射
function chaotic = tent_chaos(values)
chaotic = zeros(size(values));
chaotic(values < 0.5) = 2 * values(values < 0.5);
chaotic(values >= 0.5) = 2 * (1 - values(values >= 0.5));
end