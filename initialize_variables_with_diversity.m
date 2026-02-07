% 添加新的种群初始化方法，增加多样性
% function chromosome = initialize_variables_with_diversity(pop, M, V, min_range, max_range, item, tao, cntao, probID)
%     chromosome = initialize_variables(pop, M, V, min_range, max_range, item, tao, cntao, probID);
%     % 通过拉丁超立方体采样方法增加多样性
%     LHS_samples = lhsdesign(pop, V) .* (max_range - min_range) + min_range;
%     chromosome(:, 1:V) = LHS_samples;
% end

% % % 添加新的种群初始化方法，增加多样性（结合拉丁超立方体采样和混沌系统）
% function chromosome = initialize_variables_with_diversity(pop, M, V, min_range, max_range, item, tao, cntao, probID)
%     % 使用拉丁超立方体采样生成初始种群
%     LHS_samples = lhsdesign(pop, V) .* (max_range - min_range) + min_range;
%     chromosome = zeros(pop, M + V);
%     chromosome(:, 1:V) = LHS_samples;
% 
%     % 使用混沌系统生成一些随机种子以增加多样性
%     chaos_seeds = chaos_system(pop, V);
% 
%     % 将混沌系统生成的种子应用于初始种群
%     for i = 1:pop
%         chromosome(i, 1:V) = chromosome(i, 1:V) + chaos_seeds(i, :);
%         % 保证染色体在决策空间内
%         chromosome(i, 1:V) = max(min(chromosome(i, 1:V), max_range), min_range);
%     end
% end
% 

% 结合拉丁超立方体采样和混沌映射的初始化函数
function f = initialize_variables_with_diversity(N, M, V, min_range, max_range, item, tao, cntao, probID)
    min = min_range;
    max = max_range;
    K = M + V;

    % 使用拉丁超立方体采样生成初始种群
    LHS_samples = lhsdesign(N, V) .* (max_range - min_range) + min_range;
    
    % 使用Henon混沌映射生成混沌种子
    chaotic_samples = henon_map(N, V, 1.4, 0.3);
    chaotic_samples = chaotic_samples .* (max_range - min_range) + min_range;
    
    % 混合LHS样本和混沌样本
    f = [LHS_samples(1:N/2, :); chaotic_samples(N/2+1:end, :)];
    
    % 计算目标函数值并添加到染色体中
    for i = 1:N
        f(i,V + 1: K) = evaluate_objective(f(i,:), M, V, item, tao, cntao, probID);
    end
end

