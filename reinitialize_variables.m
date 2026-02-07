function f =reinitialize_variables(Spop,N, M, V, min_range, max_range, item, tao,cntao,radius)

%% function f = initialize_variables(N, M, V, min_tange, max_range) 
% This function initializes the chromosomes. Each chromosome has the
% following at this stage
%       * set of decision variables
%       * objective function values
% 
% where,
% N - Population size
% M - Number of objective functions
% V - Number of decision variables
% min_range - A vector of decimal values which indicate the minimum value
% for each decision variable.
% max_range - Vector of maximum possible values for decision variables.

vlb = min_range;
vub = max_range;

% K is the total number of array elements. For ease of computation decision
% variables and objective functions are concatenated to form a single
% array. For crossover(交叉) and mutation(变异) only the decision variables are used
% while for selection, only the objective variable are utilized.

K = M + V;


SPfront=paretofront(Spop,M,V);
Sfront=SPfront(:,1:V);
s=size(SPfront);

% ar_1=ceil(s(1)*rand);
% ar_2=ceil(s(1)*rand);
% direct=Sfront(ar_1,1:V)-Sfront(ar_2,1:V);
% A=null(direct,'r');     %%一列为一个方向
% A=[A direct'];


% direct=ones(1,V); 
% A=null(direct,'r');     %%一列为一个方向
% A=[A direct'];

A=eye(V);

% 
% for i=1:Hksai
%     ar=ceil(s(1)*rand);
%     d=sqrt(sum((Fpop(ar,1:V)-centroid).^2));
%      br=ceil(V*rand);
%     directnew=A(:,br)';
%      correpop(i,1:V)=Fpop(ar,1:V)+d*randn*directnew;
%     
% end

%% Initialize each chromosome
% For each chromosome perform the following (N is the population size)
for i = 1 : N
    % Initialize the decision variables based on the minimum and maximum
    % possible values. V is the number of decision variable. A random
    % number is picked between the minimum and maximum possible values for
    % the each decision variable.
   
    br=ceil(V*rand);
    directnew=A(:,br)';
%     f(i,1:V)=Spop(i,1:V)+randn*directnew*radius;
%     for j = 1 : V
% %         f(i,j) = Spop(i,j)+radius*randn;
%         if f(i,j)>vub(j)
%                f(i,j)=0.5*(Spop(i,j)+vub(j));
%         elseif f(i,j)<vlb(j)
%                 f(i,j)=0.5*(Spop(i,j)+vlb(j));
%        end
%     
%     end
ar=ceil(s(1)*rand);
    f(i,1:V)=Sfront(ar,1:V)+randn*directnew*radius;
    for j = 1 : V
%         f(i,j) = Spop(i,j)+radius*randn;
        if f(i,j)>vub(j)
               f(i,j)=0.5*(Sfront(ar,j)+vub(j));
        elseif f(i,j)<vlb(j)
                f(i,j)=0.5*(Sfront(ar,j)+vlb(j));
       end
    
    end
   

    % For ease of computation and handling data the chromosome also has the
    % vlaue of the objective function concatenated at the end. f(i,j)一行表示一个个体。The elements
    % V + 1 to K has the objective function valued. 
    % The function evaluate_objective takes one chromosome at a time,
    % infact only the decision variables are passed to the function along
    % with information about the number of objective functions which are
    % processed and returns the value for the objective functions. These
    % values are now stored at the end of the chromosome itself.
    f(i,V + 1: K) = evaluate_objective(f(i,:), M, V, item, tao,cntao);
end