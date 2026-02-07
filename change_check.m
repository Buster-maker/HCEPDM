function change = change_check(popind,M,V,item, tao,cntao,probID)
%CHANGE_CHECK Summary of this function goes here
% popind - 当前群体，一行表示一个个体，最后的M个坐标表示个体的函数值
% M - 目标数
% V - 决策变量维数
% item - 代数 tao cntao是全局变量，用以控制动态问题
s=size(popind);
index=randperm(s(1));
cum=ceil(s(1)*0.05);
count=0;
for i=1:cum
    ind=popind(index(i),:);
    %F5(x, M, V, item, tao,cntao)
    f=evaluate_objective(ind(1:V),M, V, item, tao,cntao,probID);
%        f=cec2018_DF(ind(1:V),M, V, item, tao, cntao,probID);
    diff=0;
    for j=1:M
        diff=diff+abs(ind(V+j)-f(j));
    end
    if diff<1.0E-30
        count=count+1;
    else
        break
    end
end
if count<cum
    change=true;
else
    change=false;
end
end