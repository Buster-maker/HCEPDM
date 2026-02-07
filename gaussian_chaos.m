
% 高斯混沌映射
function chaotic = gaussian_chaos(values, a)
chaotic = exp(-a .* values.^2) + 0.5;
end
