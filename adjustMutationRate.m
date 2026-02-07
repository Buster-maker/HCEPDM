% 动态调整变异率
function mutationRate = adjustMutationRate(paretoEntropy, baseRate, threshold_low, threshold_high)
    if paretoEntropy < threshold_low
        mutationRate = baseRate * 1.2; % 熵值低时增加变异率
    elseif paretoEntropy > threshold_high
        mutationRate = baseRate * 0.8; % 熵值高时减少变异率
    else
        mutationRate = baseRate; % 熵值在合理范围时保持不变
    end
end