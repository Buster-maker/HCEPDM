
% Tent映射
function chaotic = tent_chaos(values)
chaotic = zeros(size(values));
chaotic(values < 0.5) = 2 * values(values < 0.5);
chaotic(values >= 0.5) = 2 * (1 - values(values >= 0.5));
end