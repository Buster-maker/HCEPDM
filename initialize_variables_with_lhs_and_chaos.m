function f = initialize_variables_with_lhs_and_chaos(N, M, V, min_range, max_range, item, tao, cntao, probID)
% 使用拉丁超立方采样结合多级混沌映射与随机扰动的种群初始化
% N - 种群大小
% M - 目标函数数量
% V - 决策变量数量
% min_range, max_range - 决策变量范围
% item, tao, cntao, probID - 其他问题参数

% 拉丁超立方采样生成基础分布
lhs_sample = lhsdesign(N, V); % 生成[N x V]矩阵，范围在[0, 1]
min_val = min_range;
max_val = max_range;

% 将LHS样本映射到实际决策变量范围
initial_population = zeros(N, V);
for j = 1:V
    initial_population(:, j) = min_val(j) + (max_val(j) - min_val(j)) .* lhs_sample(:, j);
end

% 子群划分
num_subgroups = 3; % 划分为3个子群
subgroup_size = floor(N / num_subgroups);

% 初始化种群矩阵
f = zeros(N, M + V);

% 动态扰动比例参数
base_perturbation = 0.05; % 基础扰动比例
max_perturbation = 0.2;  % 最大扰动比例

% 遍历每个子群
for g = 1:num_subgroups
    start_idx = (g - 1) * subgroup_size + 1;
    if g == num_subgroups
        end_idx = N; % 最后一个子群包括剩余个体
    else
        end_idx = g * subgroup_size;
    end

    % 提取子群个体
    subgroup = initial_population(start_idx:end_idx, :);

    % 使用不同混沌映射扰动
    chaotic_seq = zeros(size(subgroup));
    if g == 1
        % Logistic混沌映射
        r = 3.9;
        for i = 1:size(subgroup, 1)
            chaotic_seq(i, :) = logistic_chaos(subgroup(i, :), r);
        end
    elseif g == 2
        % 高斯混沌映射
        a = 6;
        for i = 1:size(subgroup, 1)
            chaotic_seq(i, :) = gaussian_chaos(subgroup(i, :), a);
        end
    elseif g == 3
        % Tent映射
        for i = 1:size(subgroup, 1)
            chaotic_seq(i, :) = tent_chaos(subgroup(i, :));
        end
    end

    % 动态扰动比例
    dynamic_perturbation = base_perturbation + ...
        (max_perturbation - base_perturbation) * rand(size(subgroup, 1), 1);
    dynamic_perturbation = repmat(dynamic_perturbation, 1, V); % 扩展到维度V

    % 添加扰动并限制范围
    perturbed_subgroup = subgroup + dynamic_perturbation .* chaotic_seq;

    % 随机扰动增强
    random_noise = 0.01 * randn(size(perturbed_subgroup)); % 小幅随机噪声
    perturbed_subgroup = perturbed_subgroup + random_noise;

    % 保证范围合法
    perturbed_subgroup = max(min(perturbed_subgroup, max_val), min_val);

    % 保存扰动后的子群
    f(start_idx:end_idx, 1:V) = perturbed_subgroup;
end

% 计算目标函数值并附加到染色体后部
for i = 1:N
    f(i, V + 1:end) = evaluate_objective(f(i, 1:V), M, V, item, tao, cntao, probID);
end
end

% Logistic混沌映射
function chaotic = logistic_chaos(values, r)
chaotic = r .* values .* (1 - values);
end

% 高斯混沌映射
function chaotic = gaussian_chaos(values, a)
chaotic = exp(-a .* values.^2) + 0.5;
end

% Tent映射
function chaotic = tent_chaos(values)
chaotic = zeros(size(values));
chaotic(values < 0.5) = 2 * values(values < 0.5);
chaotic(values >= 0.5) = 2 * (1 - values(values >= 0.5));
end

