function preset=predictset(precenter,precenter_old,dksai,vlb,vub,M, V,Spop,item, tao,cntao,probID)

preset=[];
rowS=size(Spop);

s=dksai;
direct=sign(precenter-precenter_old);
direct_1=precenter-precenter_old;
d=sqrt(sum((precenter-precenter_old).^2));
% ar=ceil(rowS(1)*rand); 
% direct_2=Spop(ar,1:V)-precenter;
% sigama=abs(sum(direct_1.*direct_2)/(norm(direct_1)*norm(direct_2)+0.01));
for i=1:s
    ar=ceil(rowS(1)*rand);        
         preset(i,:)=Spop(ar,1:V)+direct_1+d*randn*direct;
         for j=1:V
         if preset(i,j)>vub(j)
              preset(i,j)=0.5*(vub(j)+Spop(ar,j));
        elseif  preset(i,j)<vlb(j)
                 preset(i,j)=0.5*(Spop(ar,j)+vlb(j));
         end
         end
% preset(i,:)=Spop(ar,1:V)+d*(direct_1/(norm(direct_1)+0.001)+direct_2/(norm(direct_2)+0.001))+sigama*randn*direct;
% preset(i,:)=Spop(ar,1:V)+(d+0.5*randn)*(direct_1/(norm(direct_1)+0.01)+direct_2/(norm(direct_2)+0.01));
% preset(i,:)=Spop(ar,1:V)+(direct_1+direct_2)+sigama*randn*direct;

% direct_2=Spop(ar,1:V)-precenter;
%     sigama=sum(direct_1.*direct_2)/(norm(direct_1)*norm(direct_2)+0.001);
% preset(i,:)=Spop(ar,1:V)+(direct_1+direct_2)+randn*direct; 
end

% sksai=dksai-s;
% for i=1:sksai 
%      ar=ceil(rowS(1)*rand);
%      for j=1:V
%          preset(s+i,j)=Spop(ar,j)+0.2*randn;
%          if preset(s+i,j)>vub(j)
%               preset(s+i,j)=0.5*(Spop(ar,j)+vub(j));
%         elseif  preset(s+i,j)<vlb(j)
%                  preset(s+i,j)=0.5*(Spop(ar,j)+vlb(j));
%        end
%      end
% end


for i=1:dksai
    preset(i,V + 1: V+M) = evaluate_objective(preset(i,:), M, V, item, tao,cntao,probID);
end
