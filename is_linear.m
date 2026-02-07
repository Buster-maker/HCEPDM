% 判断数据是否具有线性趋势的函数
function trend = is_linear(data)
    X = (1:length(data))';
    model = fitlm(X, data);
    R2 = model.Rsquared.Ordinary;
    threshold = 0.8; % 自定义阈值
    trend = R2 > threshold;
end