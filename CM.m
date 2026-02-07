function cm = CM(solution,M,V,T)
Aindex=find(solution(:,M+V+1)==1);
Apareto=solution(Aindex,:);
functionvalue=Apareto(:,V+1:V+M);
x=0:0.01:1;
y=1-x.^(1/((1+(0.75+0.7*sin(0.5*pi*T/10)))*(1+(0.75+0.7*sin(0.5*pi*T/10)))*10+(0.75+0.7*sin(0.5*pi*T/10)) )) ;
a=[x, y];
parameter1=reshape(a,101,2);
clear temp
popnum=size(functionvalue,1);
            truenum=size(parameter1,1);
            dis=zeros(popnum,truenum);
            fmax=max(parameter1);
            fmin=min(parameter1);
            for i=1:popnum
                for j=1:truenum
                    dis(i,j)=sqrt(sum(((functionvalue(i,:)-parameter1(j,:))./(fmax-fmin)).^2));
                end
            end
            dis=min(dis,[],2);
           cm=mean(dis);
end