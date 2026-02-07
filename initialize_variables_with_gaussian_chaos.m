function f = initialize_variables_with_gaussian_chaos(N, M, V, min_range, max_range, item, tao, cntao, probID)
% 改进后的基于高斯混沌映射的种群初始化
% 增强了种群分布的多样性和覆盖性

% 高斯混沌映射参数
a = 6; % 调整控制参数，增加混沌动态性
b = 0.5; % 偏移量，常设为0.5

min = min_range;
max = max_range;

% K 表示染色体总长度，包含决策变量和目标函数值
K = M + V;

% 初始化每个染色体
for i = 1:N
    chaotic_seq = zeros(1, V); % 用于存储混沌序列，序列长度为决策变量数V
    chaotic_seq(1) = rand(); % 混沌序列的初值随机生成

    % 使用高斯混沌映射生成混沌序列
    for j = 2:V
        chaotic_seq(j) = exp(-a * chaotic_seq(j - 1)^2) + b; % 高斯映射公式
        chaotic_seq(j) = mod(chaotic_seq(j), 1); % 确保混沌值在(0,1)范围内
    end

    % 引入随机扰动，增强分布随机性
    chaotic_seq = chaotic_seq + 0.01 * randn(1, V);
    chaotic_seq = mod(chaotic_seq, 1); % 确保扰动后仍在(0,1)范围内

    % 将混沌序列与随机序列结合
    alpha = 0.7; % 混合比例
    random_seq = rand(1, V); % 随机序列
    combined_seq = alpha * chaotic_seq + (1 - alpha) * random_seq;

    % 将混合序列映射到决策变量的定义域[min, max]
    for j = 1:V
        f(i, j) = min(j) + (max(j) - min(j)) * combined_seq(j);
    end

    % 计算目标函数值并附加到染色体后部
    f(i, V + 1:K) = evaluate_objective(f(i, :), M, V, item, tao, cntao, probID);
end
