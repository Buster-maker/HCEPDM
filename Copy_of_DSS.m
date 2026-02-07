% function DSS()  
trace1=zeros(10,50);
for run=1:10
clf;%清除当前图像窗口  
tic;%计时函数，计时的开始；toc计时的结束
pop=100;%种群中个体数目
gen=500;% 代数
tao=10;%环境变化步长、频率
item=1; %用于记录代数
ksai2=0.05; %随机迁移率
T=0;  %用于记录环境数 
ntao=10;%变化强度
cntao=1/ntao;%  
radiusold=0;%
probID='DF3';
xingneng1=[];   
baseRate = 0.05; % 基础变异率
threshold_low = 0.5; % 熵值低的阈值
threshold_high = 1.5; % 熵值高的阈值

%% 目标函数
%M是目标空间的维数，V是决策变量空间的维数。min_range和max_range是决策变量空间中变量的范围
%用户必须使用决策变量定义目标函数。objective_description_function()是另一个.m文件
    [M, V, min_range, max_range] = objective_description_function();
    centroid_old=zeros(1,V);
    center_old=zeros(1,V);
    %% Initialize the population
    %在特定范围内的随机值来进行种群的初始化，每个染色体都是由决策变量组成，同时目标函数，rank等级和拥挤距离信息的值也被添加到染色体向量中
    %但是只有当中含有函数决策变量的向量元素才会进行遗传的操作，比如交叉变异（这一步不是全部都要执行的）

    % chromosome = initialize_variables(pop, M, V, min_range, max_range, item, tao, cntao,probID);
   %高斯混沌+均匀分布
   % chromosome1 = initialize_variables_with_gaussian_chaos(pop, M, V, min_range, max_range, item, tao, cntao, probID);
chromosome = initialize_variables_with_lhs_and_chaos(pop, M, V, min_range, max_range, item, tao, cntao, probID);

    % 提取决策变量
   decision_variables= chromosome(:, 1:V);
   % decision_variables1 = chromosome1(:, 1:V);
   % decision_variables2 = chromosome2(:, 1:V);

%       % 可视化种群分布
% figure;
% 
% % 均匀随机分布种群
% subplot(1, 2, 1);
% scatter(decision_variables(:, 1), decision_variables(:, 2), 30, 'b', 'filled');
% title('均匀随机分布种群');
% xlabel('决策变量1');
% ylabel('决策变量2');
% grid on;
% xlim([min_range(1), max_range(1)]);
% ylim([min_range(2), max_range(2)]);
% 
% % % 改进的高斯混沌分布种群
% % subplot(1, 3, 2);
% % scatter(decision_variables1(:, 1), decision_variables1(:, 2), 30, 'r', 'filled');
% % title('改进的高斯混沌种群');
% % xlabel('决策变量1');
% % ylabel('决策变量2');
% % grid on;
% % xlim([min_range(1), max_range(1)]);
% % ylim([min_range(2), max_range(2)]);
% 
% % 混合混沌分布种群
% subplot(1, 2, 2);
% scatter(decision_variables2(:, 1), decision_variables2(:, 2), 30, 'g', 'filled');
% title('混合混沌种群');
% xlabel('决策变量1');
% ylabel('决策变量2');
% grid on;
% xlim([min_range(1), max_range(1)]);
% ylim([min_range(2), max_range(2)]);
% 
% % 调整图形外观
% sgtitle('两种种群初始化方法的分布对比');


    %% Sort the initialized population
    %使用非支配排序对总体进行排序。 会为每个个体返回两列，这两列的内容是所处位置对应的rank等级和拥挤度距离
    %对每条染色体来说，添加rank等级和拥挤度距离到染色体向量中，以便于计算
    chromosome = non_domination_sort_mod(chromosome, M, V);
    
    history_centroids = []; % 存储历史中心位置

    %% Start the evolution process
    % 在每一代中执行以下操作
    % * 选择适合繁殖的父代
    % * 对选定的父代执行交叉和变异的操作
    % * 从父代和其后代中进行选择
    % * 用适合的个体替换不适合的个体，以维持恒定的人口规模。
    for i = 1 : gen
        item=i;
        pool=3*pop;% pool - 交配池的大小
        tour = 2;% tour - 竞标赛大小. 二级制竞标赛选择，但是要查看这个竞标赛大小的影响，这个是任意的，是由用户决定的
%         change = change_check(chromosome,M,V,item, tao,cntao);
       mutation_rate = 0.1; 
       change = change_check(chromosome,M,V,item, tao,cntao,probID);
       
% change=0;
        if  mod(i,tao)==0%求余
         if change==true%检测环境变化
             T=T+1;

figure(1);
           switch T  %用于变量和函数的多分支选择问题
                     
               case 10
                   plot(best(:,V+1),best(:,V+2),'kO');
                   xlabel('f1'); ylabel('f2');
                   axis([0 1 0 1]);
%                            title('DSS');
                  plot(h(:,1),h(:,2),'k-');
               case 20
                    plot(best(:,V+1),best(:,V+2),'r+');
                   plot(h(:,1),h(:,2),'r-');
                  
               case 30
                   plot(best(:,V+1),best(:,V+2),'b*');
                   plot(h(:,1),h(:,2),'b-'); 
                  
               case 40
                   plot(best(:,V+1),best(:,V+2),'ms');
                   plot(h(:,1),h(:,2),'m-');
                  
               case 50
                   plot(best(:,V+1),best(:,V+2),'gp');
                   plot(h(:,1),h(:,2),'g-');
%                 
           end
           
           legend('T=10','T=10','T=20','T=20','T=30','T=30','T=40','T=40','T=50','T=50');%为图表打标注
%            legend('T=39','T=39','T=40','T=40','T=41','T=41','T=42','T=42','T=43','T=43');
           %                    legend('T=100','T=102','T=104','T=106','T=108','T=110' ,'T=112','T=114');
           hold on;

           h=PF(probID,item, tao, ntao);
           xingneng1(T)=IGD(best,M,V,item,tao,cntao,T,h);  
 trace1(run,T)= xingneng1(T);
           ParetoSF=paretofront(best,M,V); 
           ParetoS=ParetoSF(:,1:V);
           [row ,temp]=size(ParetoS);
%            clear temp
           center_new=sum(ParetoS,1)./row;


history_centroids = [history_centroids; center_new];

% 使用最适合的模型预测下一代的中心位置
% 初始化预测的下一代中心位置
next_center = zeros(1, V);
weights = zeros(1, 2); % 存储ARIMA和SVR的权重

% 对每个决策变量维度进行预测
for d = 1:V
    % 如果历史中心位置记录数量大于2
    if size(history_centroids, 1) > 2
        % 选择最适合的模型（线性或非线性）
        trend = is_linear(history_centroids(:, d));
        if trend
            % 如果是线性趋势，使用ARIMA模型进行预测
            try
                model_arima = arima('Constant',0,'D',1,'Seasonality',0);
                fruoarima = estimate(model_arima, history_centroids(:, d), 'Display', 'off');
                [forecast_arima, ~] = forecast(fit_arima, 1);
                next_center(d) = forecast_arima;
                weights = [1, 0]; % 完全依赖ARIMA
            catch
                next_center(d) = history_centroids(end, d);
            end
        else
            % 如果是非线性趋势，使用ARIMA和SVR模型进行预测并动态调整权重
            try
                % ARIMA模型预测
                model_arima = arima('Constant',0,'D',1,'Seasonality',0);
                fit_arima = estimate(model_arima, history_centroids(:, d), 'Display', 'off');
                [forecast_arima, ~] = forecast(fit_arima, 1);
            catch
                forecast_arima = history_centroids(end, d);
            end

            try
                % 支持向量机回归（SVR）预测
                model_svr = fitrsvm((1:size(history_centroids, 1))', history_centroids(:, d), 'KernelFunction', 'gaussian', 'KernelScale', 'auto');
                forecast_svr = predict(model_svr, size(history_centroids, 1) + 1);
            catch
                forecast_svr = history_centroids(end, d);
            end

            % 计算最近几次预测的误差
            error_arima = abs(forecast_arima - history_centroids(end, d));
            error_svr = abs(forecast_svr - history_centroids(end, d));
            total_error = error_arima + error_svr;
            if total_error > 0
                weights(1) = error_svr / total_error; % 权重与另一模型的误差成正比
                weights(2) = error_arima / total_error;
            else
                weights = [0.5, 0.5]; % 如果总误差为零，则平均权重
            end

            % 动态加权平均
            next_center(d) = weights(1) * forecast_arima + weights(2) * forecast_svr;
        end
    else
        % 如果历史中心位置记录数量小于等于2，直接使用当前中心位置
        next_center(d) = center_new(d);
    end
end



            %%更新种群位置
             movingdirect= center_new- center_old;  
             radiusnew=norm(movingdirect);

             radius=max( [radiusnew,radiusold] );
             radiusold=radiusnew;
           
            chromosome = reinitialize_predict(chromosome,pop, M, V, min_range, max_range, item, tao,cntao,movingdirect,probID);%DSS1
            % chromosome = reinitialize_predict_gaussuian(chromosome,pop, M, V, min_range, max_range, item, tao,cntao,movingdirect,probID);%DSS1
          
     
            chromosome = non_domination_sort_mod(chromosome, M, V);
           centroid_old=zeros(1,V);   
    %        end
           center_old=center_new;
           
         end

        end
    % 选择过程
    %二进制随机数和它们的适应度进行比较，适应度高的个体被选为作为父代，比赛选择一直到交配池满为止，
    %一般来说池的大小就是被选择作为父代的个数
    %函数tournament_selection输入的参数是 染色体、交配池大小、竞标赛大小
        parent_chromosome = tournament_selection(chromosome, pool, tour);


        % 进行交叉和变异操作
        %交叉算子和变异算子的分布指数分别为mu = 20和mum = 20
    %     mu = 10;    
        mum = 20;
        offspring_chromosome = ...
            genetic_operator(parent_chromosome, ...
            M, V, mum, min_range, max_range, item, tao,cntao,probID);

        % Intermediate population（中间种群）
    %中间种群是当前这一代父母和后代的总人口。 人口规模是初始人口的两倍。

        [main_pop,temp] = size(chromosome);
        [offspring_pop,temp] = size(offspring_chromosome);
        intermediate_chromosome(1:main_pop,:) = chromosome;
        intermediate_chromosome(main_pop + 1 : main_pop + offspring_pop,1 : M+V) = ...
            offspring_chromosome;

        % Non-domination-sort of intermediate population中间种群的非支配排序
        intermediate_chromosome = ...
            non_domination_sort_mod(intermediate_chromosome, M, V);
        % Perform Selection执行选择
        chromosome = replace_chromosome(intermediate_chromosome, M, V, pop);   


        best=chromosome;
       ParetoSF=paretofront(chromosome,M,V); 
        ParetoS=ParetoSF(:,1:V);

% 
% 计算 Pareto 熵
 [paretoEntropy, deltaEntropy] = calculateParetoEntropy(ParetoS);

 % 动态调整变异率
mutationRate = adjustMutationRate(paretoEntropy, baseRate, threshold_low, threshold_high);

% 变异操作以增加种群多样性
chromosome = mutatePopulation(chromosome, mutationRate, min_range, max_range, paretoEntropy, probID);

        dksai=floor(ksai2*pop);
        [row ,temp]=size(ParetoS);
        clear temp
        centroid=sum(ParetoS,1)./row;
        preset=predictset(centroid,centroid_old,dksai,min_range,max_range,M, V,ParetoS,item, tao,cntao,probID);
        centroid_old=[];
        centroid_old=centroid; 
        suijixuanqu=randperm(pop);
        chromosome(suijixuanqu(1:dksai),1:M+V)=preset;
        chromosome=non_domination_sort_mod(chromosome, M, V);
if ~mod(item,10)
        disp([num2str(item),' generations completed, chromesome number ',num2str(size(chromosome,1))])
end  
  
toc
 end    
 if ~mod(run,2)
        disp([num2str(run),' run completed, chromesome number ',num2str(size(run,1))])
 end   
 end  
%%以下完全原创
for i=1:50
    A(i)=mean(trace1(:,i));
end

figure(2); 
plot(A,'b-');
%  saveas(H, 'filename','format')
% H表示图像的句柄，format的选择有fig，m，mfig，mmat，jpg 等
saveas(figure(1),num2str(probID),'fig');
saveas(figure(1),num2str(probID),'tif');
saveas(figure(2),num2str(probID),'fig');

MIGD=mean(mean(trace1));
MIGD_1=mean(mean(trace1(:,1:10)));
MIGD_2=mean(mean(trace1(:,11:40)));
MIGD_3=mean(mean(trace1(:,41:50)));
TRACE=trace1.';
SD=mean(std(TRACE));
SD_1=mean(std(TRACE(1:10,:)));
SD_2=mean(std(TRACE(11:40,:)));
SD_3=mean(std(TRACE(41:50,:)));
MIGD=mean(mean(trace1))
MIGD_1=mean(mean(trace1(:,1:10)))
MIGD_2=mean(mean(trace1(:,11:40)))
MIGD_3=mean(mean(trace1(:,41:50)))
SD=mean(std(TRACE))
SD_1=mean(std(TRACE(1:10,:)))
SD_2=mean(std(TRACE(11:40,:)))
SD_3=mean(std(TRACE(41:50,:)))
save([(regexprep(probID, '[^\w'']', '_')) '.mat']);
%% Result
% Save the result in ASCII text format.
% save solution.txt best -ASCII

%% Visualize
% The following is used to visualize the result if objective space
% dimension is visualizable.
% if M == 2
%     plot(chromosome(:,V + 1),chromosome(:,V + 2),'*');
% elseif M ==3
%     plot3(chromosome(:,V + 1),chromosome(:,V + 2),chromosome(:,V + 3),'*');
% end
  
