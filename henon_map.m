% Henon混沌映射生成器
function chaotic_sequence = henon_map(N, V, a, b)
    chaotic_sequence = zeros(N, V);
    % 初始化x和y
    x = rand(N, 1);
    y = rand(N, 1);
    for i = 2:N
        x(i) = 1 - a * x(i-1)^2 + y(i-1);
        y(i) = b * x(i-1);
    end
    chaotic_sequence(:,1) = x;
    chaotic_sequence(:,2:end) = repmat(x, 1, V-1) + 0.01 * rand(N, V-1);
end