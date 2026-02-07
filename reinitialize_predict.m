function f =reinitialize_predict(Spop,N, M, V, min_range, max_range, item, tao,cntao,movingdirect,probID)

%% function f = initialize_variables(N, M, V, min_tange, max_range) 

vlb = min_range;
vub = max_range;

K = M + V;
SPfront=paretofront(Spop,M,V); 
Sfront=SPfront(:,1:V);
s=size(SPfront);

direct=movingdirect;
d=norm(direct);
A=null(direct,'r');     %%一列为一个方向
signdirect=sign(direct);
N1=N/2;
N2=N-N1;
for i = 1 : N1
ar=ceil(s(1)*rand);
f(i,1:V)=Sfront(ar,1:V)+direct+d*randn*signdirect;
    for j = 1 : V
        if f(i,j)>vub(j)
               f(i,j)=0.5*(Sfront(ar,j)+vub(j));
        elseif f(i,j)<vlb(j)
                f(i,j)=0.5*(Sfront(ar,j)+vlb(j));
       end
    
    end

    f(i,V + 1: K) = evaluate_objective(f(i,:), M, V, item, tao,cntao,probID);
end

%% Initialize each chromosome
for i = 1 : N2
  
    br=ceil((V-1)*rand);
    directnew=A(:,br)';
    
    
%     f(i,1:V)=Spop(i,1:V)+directnew+randn*signdirect;
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
f(N1+i,1:V)=Sfront(ar,1:V)+randn*directnew;
    for j = 1 : V
%         f(i,j) = Spop(i,j)+radius*randn;
        if f(N1+i,j)>vub(j)
               f(N1+i,j)=0.5*(Sfront(ar,j)+vub(j));
        elseif f(N1+i,j)<vlb(j)
                f(N1+i,j)=0.5*(Sfront(ar,j)+vlb(j));
       end
    
    end
    f(N1+i,V + 1: K) = evaluate_objective(f(N1+i,:), M, V, item, tao,cntao,probID);
end