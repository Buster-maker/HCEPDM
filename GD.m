function gd = GD(solution,M,V,T)
%IGD Summary of this function goes here
%   Detailed explanation goes here

Aindex=find(solution(:,M+V+1)==1);
Apareto=solution(Aindex,:);
[r,temp]=size(Apareto);
clear temp
x=0:0.01:1;
y=1-x.^(1/((1+(0.75+0.7*sin(0.5*pi*T/10)))*(1+(0.75+0.7*sin(0.5*pi*T/10)))*10+(0.75+0.7*sin(0.5*pi*T/10)) )) ;
dmin=[];
for i=1:r
    b=Apareto(i,V+1:V+M);
    for j=1:100
        a=[x(j), y(j)];
        d(j)=sqrt(sum((b-a).^2));
    end
    dmin(i)=min(d);
    d=[];
end
gd=sum(dmin)/r;



    





    
