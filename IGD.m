function igd = IGD(chromosome,M,V,item,tao,cntao,T,h)
%IGD Summary of this function goes here
%   Detailed explanation goes here

Aindex=find(chromosome(:,M+V+1)==1);
Apareto=chromosome(Aindex,:);
[r,temp]=size(Apareto);
clear temp

[test,~]=size(h);
% x=0:0.002:1;
% H=1.25+0.75*sin(0.5*pi*cntao*(T-1));
% y=1-x.^(1/((1+(0.75+0.7*sin(0.5*pi*T/10)))*(1+(0.75+0.7*sin(0.5*pi*T/10)))*10+(0.75+0.7*sin(0.5*pi*T/10)) )) ;
% y=1-x.^(1/((1+(1.25+0.75*sin(0.5*pi*T/10)))*(1+(1.25+0.75*sin(0.5*pi*T/10)))*10+(1.25+0.7*sin(0.5*pi*T/10)) )) ;
% y=1-sqrt(x);%fda1

% x=s.^H;
% y=(1-s).^H;
% x=s;
% y=1-x.^H;
dmin=[];
for i=1:test
    b=[h(i,1), h(i,2)];
    for j=1:r
        a=Apareto(j,V+1:V+M);
        d(j)=sqrt(sum((b-a).^2));
    end
    dmin(i)=min(d);
    d=[];
end
igd=sum(dmin)/test;
    
end

% function igd = IGD( solution,M,V,T)
% %IGD Summary of this function goes here
% %   Detailed explanation goes here
% 
% Aindex=find(solution(:,M+V+1)==1);
% Apareto=solution(Aindex,:);
% [r,temp]=size(Apareto);
% clear temp
% x=0:0.01:1;
% y=1-x.^(1/((1+(0.75+0.7*sin(0.5*pi*T/10)))*(1+(0.75+0.7*sin(0.5*pi*T/10)))*10+(0.75+0.7*sin(0.5*pi*T/10)) )) ;
% % y=x.^(1/2);
% dmin=[];
% for i=1:101
%     b=[x(i), y(i)];
%     for j=1:r
%         a=Apareto(j,V+1:V+M);
%         d(j)=sqrt(sum((b-a).^2));
%     end
%     dmin(i)=min(d);
%     d=[];
% end
% igd=sum(dmin)/101;

% y=1-x.^(1/(H*(1+H)*10+H)) ;
% y=1-x.^(1/((1+(0.75+0.7*sin(0.5*pi*T/10)))*(1+(0.75+0.7*sin(0.5*pi*T/10)))*10+(0.75+0.7*sin(0.5*pi*T/10)) )) ;



