function f = reinitialize_predict_gaussuian(Spop, N, M, V, min_range, max_range, item, tao, cntao, movingdirect, probID)

% 设置比例
ratioChaos = 0.6; % 高斯混沌比例
ratioRandom = 0.4; % 均匀分布比例

popChaos = round(N * ratioChaos);
popRandom = N - popChaos;

% 高斯混沌映射初始化
a = 1.4;
b = 0.3;
x = rand;
y = rand;
chaotic_sequence_x = zeros(popChaos * V, 1);
chaotic_sequence_y = zeros(popChaos * V, 1);
chaotic_sequence_x(1) = x;
chaotic_sequence_y(1) = y;
for i = 2:popChaos * V
    chaotic_sequence_x(i) = 1 - a * chaotic_sequence_x(i-1)^2 + chaotic_sequence_y(i-1);
    chaotic_sequence_y(i) = b * chaotic_sequence_x(i-1);
end
chaos_samples = zeros(popChaos, V);
for i = 1:popChaos
    for j = 1:V
        index = (i - 1) * V + j;
        chaos_samples(i, j) = min_range(j) + chaotic_sequence_x(index) * (max_range(j) - min_range(j));
        chaos_samples(i, j) = max(min(chaos_samples(i, j), max_range(j)), min_range(j));
    end
end

% 均匀分布初始化
random_samples = zeros(popRandom, V);
for i = 1:popRandom
    for j = 1:V
        random_samples(i, j) = min_range(j) + (max_range(j) - min_range(j)) * rand();
    end
end

% 合并所有样本
decision_variables = [chaos_samples; random_samples];

% 初始化种群
SPfront = paretofront(Spop, M, V); 
Sfront = SPfront(:, 1:V);
s = size(SPfront);

direct = movingdirect;
d = norm(direct);
A = null(direct, 'r');
signdirect = sign(direct);
N1 = N / 2;
N2 = N - N1;

% 初始化前 N1 个个体
for i = 1 : N1
    ar = ceil(s(1) * rand);
    f(i, 1:V) = Sfront(ar, 1:V) + direct + d * randn * signdirect;
    for j = 1 : V
        if f(i, j) > max_range(j)
            f(i, j) = 0.5 * (Sfront(ar, j) + max_range(j));
        elseif f(i, j) < min_range(j)
            f(i, j) = 0.5 * (Sfront(ar, j) + min_range(j));
        end
    end
    f(i, V + 1:V + M) = evaluate_objective(f(i, 1:V), M, V, item, tao, cntao, probID);
end

% 使用高斯混沌映射和均匀分布重新初始化后 N2 个个体
for i = 1:N2
    if i <= popChaos
        f(N1 + i, 1:V) = chaos_samples(i, :);
    else
        f(N1 + i, 1:V) = random_samples(i - popChaos, :);
    end
    for j = 1:V
        if f(N1 + i, j) > max_range(j)
            f(N1 + i, j) = 0.5 * (f(N1 + i, j) + max_range(j));
        elseif f(N1 + i, j) < min_range(j)
            f(N1 + i, j) = 0.5 * (f(N1 + i, j) + min_range(j));
        end
    end
    f(N1 + i, V + 1:V + M) = evaluate_objective(f(N1 + i, 1:V), M, V, item, tao, cntao, probID);
end

end
