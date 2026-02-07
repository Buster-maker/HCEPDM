% 变异操作函数
% function chromosome = mutatePopulation(chromosome, mutationRate, min_range, max_range)
%     [pop, K] = size(chromosome);
%     V = length(min_range); % 决策变量数
%     for i = 1:pop
%         if rand < mutationRate
%             % 对每个个体进行变异
%             for j = 1:V
%                 % 随机选择一个基因进行变异
%                 if rand < mutationRate
%                     chromosome(i, j) = min_range(j) + (max_range(j) - min_range(j)) * rand();
%                 end
%             end
%         end
%     end
% end
% function chromosome = mutatePopulation(chromosome, mutationRate, min_range, max_range, paretoEntropy)
%     [pop, K] = size(chromosome);
%     V = length(min_range); % 决策变量数
%     for i = 1:pop
%         if rand < mutationRate
%             % 随机选择一个基因进行变异
%             for j = 1:V
%                 if rand < mutationRate * paretoEntropy % 根据 Pareto 熵调整基因变异率
% 
%                     %chromosome(i, j) = min_range(j) + (max_range(j) - min_range(j)) * rand();
% 
%                     % current_value = chromosome(i, j);
%                     % %将当前基因值归一化到[0, 1]范围内
%                     % normalized_value = (current_value - min_range(j)) / (max_range(j) - min_range(j));
%                     % %将归一化后的值重新映射到[min_range, max_range]范围内
%                     % chromosome(i, j) = min_range(j) + (max_range(j) - min_range(j)) * normalized_value;
% 
%                 end
%             end
%         end
%     end
% end

function chromosome = mutatePopulation(chromosome, mutationRate, min_range, max_range, paretoEntropy, probID)
%     [pop, K] = size(chromosome);
%     V = length(min_range); % 决策变量数
% 
%     % 获取需要调整的维度和其他参数
%     CAJ = CaJ_erwei(probID);
% 
%     for i = 1:pop
%         if rand < mutationRate
%             % 随机选择一个基因进行变异
%             for j = 1:V
%                 if rand < mutationRate * paretoEntropy % 根据 Pareto 熵调整基因变异率
%                     % 将当前基因值进行均匀分布调整
%                     if j == CAJ(1, 2)
%                         max1 = max(chromosome(:, j));
%                         min1 = min(chromosome(:, j));
%                         % 将当前基因值调整为均匀分布
%                         BY = min1 + ((max1 - min1) / pop) * i;
%                         chromosome(i, j) = max(min(BY, max_range(j)), min_range(j));
%                     else
%                         % 对其他基因进行标准的随机变异
%                         current_value = chromosome(i, j);
%                         normalized_value = (current_value - min_range(j)) / (max_range(j) - min_range(j));
%                         chromosome(i, j) = min_range(j) + (max_range(j) - min_range(j)) * normalized_value;
%                     end
%                 end
%             end
%         end
%     end
% end

    [pop, K] = size(chromosome);
    V = length(min_range); % 决策变量数

    % 获取需要调整的维度和其他参数
    CAJ = CaJ_erwei(probID);

    % 混沌映射参数
    r = 3.9; % Logistic映射参数
    a = 6;   % 高斯映射参数
    tent_threshold = 0.5; % Tent映射分界点

    for i = 1:pop
        if rand < mutationRate
            % 随机选择一个基因进行变异
            for j = 1:V
                if rand < mutationRate * paretoEntropy % 根据 Pareto 熵调整基因变异率
                    % 混沌扰动
                    chaotic_value = 0;
                    chaotic_choice = randi([1, 3]); % 随机选择一种混沌映射
                    if chaotic_choice == 1
                        % Logistic混沌映射
                        chaotic_value = r * rand * (1 - rand);
                    elseif chaotic_choice == 2
                        % 高斯混沌映射
                        chaotic_value = exp(-a * rand^2) + 0.5;
                    elseif chaotic_choice == 3
                        % Tent混沌映射
                        chaotic_value = rand;
                        if chaotic_value < tent_threshold
                            chaotic_value = 2 * chaotic_value;
                        else
                            chaotic_value = 2 * (1 - chaotic_value);
                        end
                    end

                    % 应用混沌扰动
                    if j == CAJ(1, 2)
                        % 特定维度使用混沌值扰动
                        max1 = max(chromosome(:, j));
                        min1 = min(chromosome(:, j));
                        BY = min1 + chaotic_value * (max1 - min1);
                        chromosome(i, j) = max(min(BY, max_range(j)), min_range(j));
                        
                    else
                        % 其他维度的随机变异结合混沌扰动
                        current_value = chromosome(i, j);
                        normalized_value = (current_value - min_range(j)) / (max_range(j) - min_range(j));
                        chaotic_adjustment = normalized_value + chaotic_value * (1 - normalized_value);
                        chromosome(i, j) = min_range(j) + (max_range(j) - min_range(j)) * chaotic_adjustment;
                    end
                end
            end
        end
    end
end


