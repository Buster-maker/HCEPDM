% Logistic混沌映射
function chaotic = logistic_chaos(values, r)
chaotic = r .* values .* (1 - values);
end
