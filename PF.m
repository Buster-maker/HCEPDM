function h=PF(probID,item, tao, ntao)
% T0=tao;
% calculate time instant
% tau_tmp=max(item+item-(T0+1),0);
t=1/ntao*floor(item/tao);
g=1;
% H=20; % number of divisions along each objective.
V=20;
switch (probID)
    case 'FDA1'
        x=linspace(0,1,1500);
        f1=x;
        H=sin(0.5*pi*t);
        f2=1-(x/g).^(1/2);
        [h]=get_PF({f1,f2}, false);
    case 'FDA2'
        x=linspace(0,1,1500);
        f1=x;
        H=0.75+0.7*sin(0.5*pi*t);
        x3s=sum((x(V/2+1:V)-H).^2);
        z=1-(x/g).^(1/(H+x3s));
        f2=g*z;
       [h]=get_PF({f1,f2}, false);
%     case 'FDA3'
%         x=linspace(0,1,1500);
%         F=10^(2*sin(0.5*pi*t));
%         f1=x.^F;
%         G=abs(sin(0.5*pi*t));
%         z=1-(x/g).^(1/2);
%         f2=g*z;
%         [h]=get_PF({f1,f2}, false);
case 'FDA3'
        x=linspace(0,1,1000);
        F=10^(2*sin(0.5*pi*t));        
        G=abs(sin(0.5*pi*t));
        f1=(x.^F);
        f2=(1+G)*(1-sqrt(f1./(1+G)));        
        [h]=get_PF({f1,f2}, false);
    case 'FDA4'
    case 'DMOP1'
        x=linspace(0,1,1500);
        f1=x;
        H=1.25+0.75*sin(0.5*pi*t);
        z=1-(x/g).^(H);
        f2=g*z;
        [h]=get_PF({f1,f2}, false);
        
        
     case 'DMOP2'   
       x=linspace(0,1,1500);
        H=0.75*sin(0.5*pi*t)+1.25;
        f1=x;
        f2=g*(1-(x/g).^H);
        [h]=get_PF({f1,f2}, false); 
		
	case 'DMOP3'     
        x=linspace(0,1,1500);
        f1=x;
        f2=g*(1-(x/g).^0.5);
        [h]=get_PF({f1,f2}, false); 
      
    case 'F5'
        x=linspace(0,1,1500);
        %         f=[];y=[];
        H=1.25+0.75*sin(pi*t);
        %         fts=t-floor(t);
        
        %         a=2*(cos(fts*pi*0.1)+1);
        %         b=2*(sin(2*t*pi)+1);
        %         for i=1:V
        %             y(i)=x(i)-b-1+(abs(x(1)-a))^(H+i/V);
        %         end
        %         f1=x.^H+1;
        %         f2=(1-x).^H+1;
        
        f1=(x).^H;
        f2=(1-x).^H;
        [h]=get_PF({f1,f2}, false);
    case 'F6'
        x=linspace(0,1,1500);
        %         f=[];y=[];
        H=1.25+0.75*sin(pi*t);
        %         fts=t-floor(t);
        %         fts=t-floor(t);
        %         a=2*(cos(fts*pi*0.1)+1);
        %         b=2*(sin(2*t*pi)+1);
        %         for i=1:V
        %             y(i)=x(i)-b-1+(abs(x(1)-a))^(H+i/V);
        %         end
        %         f1=x.^H+1;
        %         f2=(1-x).^H+1;
        
        f1=(x).^H;
        f2=(1-x).^H;
        [h]=get_PF({f1,f2}, false);
            case 'F7'
        x=linspace(0,1,1500);
        %         f=[];y=[];
        H=1.25+0.75*sin(pi*t);
        %         fts=t-floor(t);
        %         fts=t-floor(t);
        %         a=2*(cos(fts*pi*0.1)+1);
        %         b=2*(sin(2*t*pi)+1);
        %         for i=1:V
        %             y(i)=x(i)-b-1+(abs(x(1)-a))^(H+i/V);
        %         end
        %         f1=x.^H+1;
        %         f2=(1-x).^H+1;
        
        f1=(x).^H;
        f2=(1-x).^H;
        [h]=get_PF({f1,f2}, false);
    case 'DF1'   %dmop2
        x=linspace(0,1,1500);
        H=0.75*sin(0.5*pi*t)+1.25;
        f1=x;
        f2=g*(1-(x/g).^H);
        [h]=get_PF({f1,f2}, false);
    case 'DF2'   %¸Ä½ødmop3
        x=linspace(0,1,1500);
        G=abs(sin(0.5*pi*t));
        f1=x;
        f2=g*(1-(x/g).^0.5);
        [h]=get_PF({f1,f2}, false);
    case 'DF3'   %zjz
        x=linspace(0,1,1500);
        G=sin(0.5*pi*t);
        H=G+1.5;
        f1=x;
        f2=g*(1-(x/g).^H);
        [h]=get_PF({f1,f2}, false);
    case 'DF4'
        a=sin(0.5*pi*t);
        b=1+abs(cos(0.5*pi*t));     
        x=linspace(a,a+b);
        H=1.5+a;
        f1=g*abs(x-a).^H;
        f2=g*abs(x-a-b).^H;
        [h]=get_PF({f1,f2}, false);
    case 'DF5'
        x=linspace(0,1,1500);
        G=sin(0.5*pi*t);
        w=floor(10*G);
        f1=g*(x+0.02*sin(w*pi*x));
        f2=g*(1-x+0.02*sin(w*pi*x));
        [h]=get_PF({f1,f2}, false);
    case 'DF6'
        x=linspace(0,1,1500);
        G=sin(0.5*pi*t);
        a=0.2+2.8*abs(G);
        f1=g*(x+0.1*sin(3*pi*x)).^a;
        f2=g*(1-x+0.1*sin(3*pi*x)).^a;
        [h]=get_PF({f1,f2}, false);
    case 'DF7'
        x=linspace(1,4,1500);
        f1=g*(1+t)./x;
        f2=g*x/(1+t);
        [h]=get_PF({f1,f2}, false);
    case 'DF8'
        x=linspace(0,1,1500);
        a=2.25+2*cos(2*pi*t);
        f1=g*(x+0.1*sin(3*pi*x));
        f2=g*(1-x+0.1*sin(3*pi*x)).^a;
        [h]=get_PF({f1,f2}, false);
    case 'DF9'
        x=linspace(0,1,1500);
        N=1+floor(10*abs(sin(0.5*pi*t)));
        f1=g*(x+max(0, (0.1+0.5/N)*sin(2*N*pi*x)));
        f2=g*(1-x+max(0, (0.1+0.5/N)*sin(2*N*pi*x)));
        [h]=get_PF({f1,f2}, true); 
    case 'DF10'
        [x1,x2]=meshgrid(linspace(0,1,H));
        H=2.25+2*cos(0.5*pi*t);
        f1=g*sin(0.5*pi*x1).^H;
        f2=g*sin(0.5*pi*x2).^H.*cos(0.5*pi*x1).^H;
        f3=g*cos(0.5*pi*x2).^H.*cos(0.5*pi*x1).^H;
        [h]=get_PF({f1,f2, f3}, false); % PF is continous, so no need for nondominated sorting
    case 'DF11'
        [x1,x2]=meshgrid(linspace(0,1,H));
        G=abs(sin(0.5*pi*t));
        y1=pi*G/6+(pi/2-pi*G/3)*x1;
        y2=pi*G/6+(pi/2-pi*G/3)*x2;
        f1=g.*sin(y1);
        f2=g.*sin(y2)*cos(y1);
        f3=g.*cos(y2)*cos(y1);
        [h]=get_PF({f1,f2, f3}, false); % PF is continous, so no need for nondominated sorting
    case 'DF12'
        [x1,x2]=meshgrid(linspace(0,1,H));
        k=10*sin(pi*t);
        tmp2=abs(sin(floor(k*(2*x1-1))*pi/2).*sin(floor(k*(2*x2-1))*pi/2));
        g=1+tmp2;
        f1=g.*cos(0.5*pi*x2).*cos(0.5*pi*x1);
        f2=g.*sin(0.5*pi*x2).*cos(0.5*pi*x1);
        f3=g.*sin(0.5*pi*x1);
        [h]=get_PF({f1,f2, f3}, true); 
    case 'DF13'
        [x1,x2]=meshgrid(linspace(0,1,H));
        G=sin(0.5*pi*t);
        p=floor(6*G);
        f1=g.*cos(0.5*pi*x1).^2;
        f2=g.*cos(0.5*pi*x2).^2;
        f3=g.*sin(0.5*pi*x1).^2+sin(0.5*pi*x1).*cos(p*pi*x1).^2 +...
            sin(0.5*pi*x2).^2+sin(0.5*pi*x2).*cos(p*pi*x2).^2;
        [h]=get_PF({f1,f2, f3}, true);
    case 'DF14'
        [x1,x2]=meshgrid(linspace(0,1,H));
        G=sin(0.5*pi*t);
        y=0.5+G*(x1-0.5);
        f1=g.*(1-y+0.05*sin(6*pi*y));
        f2=g.*(1-x2+0.05*sin(6*pi*x2)).*(y+0.05*sin(6*pi*y));
        f3=g.*(x2+0.05*sin(6*pi*x2)).*(y+0.05*sin(6*pi*y));
        [h]=get_PF({f1,f2, f3}, false); % PF is continous, so no need for nondominated sorting
    otherwise
        disp('no such test problem.')
end
end
% helper function: identify nondominated solutions
% f={f1,f2,f3,...}
% nondominate: nondominated sorting (true), otherwise (false)
function [h]=get_PF(f, nondominate)
ncell=length(f);
s=numel(f{1});
h=[];
for i=1:ncell
    fi=reshape(f{i},s,1);
    h=[h,fi];
end

if nondominate
    in=get_skyline(h);
    h=h(in,:);
end
end


% helper function: find the indices of nondominated solutions
function [Im]=get_skyline(x)

n=size(x,1);

if n==1
    Im=1;
else
    Iold=1:n;
    Im=[]; %Im holds the index of the solutions in the first front;
    
    % sieving method
    %[x1_index,x1]=min(sum(x'));
    z=x;
    n1=n-1;
    while (n1>0)
        [V]=get_vector(z);
        I=[];
        J=[];
        %J=find(V==1);
        %I=find(x<0,1);
        for k=1:n1
            if V(k)==1
                J=[J,k];
            elseif V(k)==-1
                I=[I,k];
            end
        end
        
        if isempty(I)
            Im=[Im,Iold(1)];
           % plot(x(Iold(1),1),x(Iold(1),2),'ro');
           % hold on;
        end
        In=setdiff(2:n1+1,J+1);
        
        In=Iold(In);
        if length(In)==1
            Im=[Im,In];
            break;
        end
        
        z=x(In,:);
        n1=length(In)-1; %because the size of get_vector(z)is smaller than In.
        Iold=In;
    end
end
end

function [myvector,ob_count]=get_vector(x)
ob_count=0;
[N,D]=size(x);

myvector=[];
for j=2:N
    dom_less=0;
    dom_equal=0;
    dom_more=0;
    %p=sign( x(1,:)-x(2,:))
    for k = 1 : D
        if (x(1,k) < x(j,k))
            dom_less = dom_less + 1;
        elseif (x(1, k) == x(j,k))
            dom_equal = dom_equal + 1;
        else
            dom_more = dom_more + 1;
        end
        ob_count=ob_count+1;
        if dom_less>1&&dom_more>1
            break;
        end
    end
    if k<D
        myvector(j-1)=0;
        continue;
    end
    if dom_less == 0 && dom_equal ~= D
        myvector(j-1)=-1;
    elseif dom_more == 0 && dom_equal ~= D
        myvector(j-1)=1;
    else
        myvector(j-1)=0;
    end
end
end


