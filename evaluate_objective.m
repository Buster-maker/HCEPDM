function f = evaluate_objective(x, M, V, item, tao, cntao, probID)
t=cntao*floor(item/tao);
f=[];
n=V;
switch (probID)
    case 'FDA1'
        f(1)=x(1);
        H=sin(0.5*pi*t);
        g=1+sum((x(2:V)-H).^2);
        h=1-(f(1)/g).^(1/2);
        f(2)=g*h;
    case 'FDA2'
        f(1)=x(1);
        H=0.75+0.7*sin(0.5*pi*t);
        g=1+sum(x(2:V/2).^2);
        x3s=sum((x(V/2+1:V)-H).^2);
        h=1-(f(1)/g).^(1/(H+x3s));
        f(2)=g*h;
    case 'FDA3'
        F=10^(2*sin(0.5*pi*t));
        f(1) =x(1).^F;
        G=abs(sin(0.5*pi*t));
        g=1+G+sum((x(2:V)-G).^2);
        h=1-(f(1)/g).^(1/2);
        f(2)=g*h;
    case 'FDA4'
    case 'DMOP1'
        f(1)=x(1);
        g=1+9*sum(x(2:V).^2);
        H=1.25+0.75*sin(0.5*pi*t);
        h=1-(f(1)/g).^(H);
        f(2)=g*h;
   case 'DMOP2'     
         f(1)=x(1);
        G=abs(sin(0.5*pi*t));
        H=1.25+0.75*sin(0.5*pi*t);
        g=1+sum((x(2:V)-G).^2);
        f(2)=g*(1-(f(1)/g).^H);
	case 'DMOP3'   
        f(1)=x(1);
        G=abs(sin(0.5*pi*t));
        g=1+sum((x(2:V)-G).^2);
        f(2)=g*(1-(f(1)/g).^0.5);
    case 'F5'
        f=[];y=[];                
        H=1.25+0.75*sin(pi*t*0.1);
%         fts=t-floor(t);
        a=2*(cos(t*pi)+1);
        b=2*(sin(2*t*pi)+1);
        for i=1:V
            y(i)=x(i)-b-1+(abs(x(1)-a))^(H+i/V);
        end
        f(1)=(abs(x(1)-a))^H+sum(y(3:2:V-1).^2);
        f(2)=(abs(x(1)-a-1))^H+sum(y(2:2:V).^2);
    case 'F6'
         f=[];y=[];                
        H=1.25+0.75*sin(pi*t);
%         fts=t-floor(t);
        a=2*cos(1.5*t*pi)*sin(0.5*pi*t)+2;
        b=2*(sin(1.5*t*pi)*cos(0.5*pi*t)+1);
        for i=1:V
            y(i)=x(i)-b-1+(abs(x(1)-a))^(H+i/V);
        end
        f(1)=(abs(x(1)-a))^H+sum(y(3:2:V-1).^2);
        f(2)=(abs(x(1)-a-1))^H+sum(y(2:2:V).^2);
    case 'F7'
        f=[];y=[];
        H=1.25+0.75*sin(pi*t);
        %         fts=t-floor(t);
        a=1.7*(1-sin(pi*t))*sin(pi*t)+3.4;
        b=1.4*(1-sin(t*pi))*cos(pi*t)+2.1;
        for i=1:V
            y(i)=x(i)-b-1+(abs(x(1)-a))^(H+i/V);
        end
        f(1)=(abs(x(1)-a))^H+sum(y(3:2:V-1).^2);
        f(2)=(abs(x(1)-a-1))^H+sum(y(2:2:V).^2);
    case 'DF1'
        %DMOP2
        f(1)=x(1);
        G=abs(sin(0.5*pi*t));
        H=1.25+0.75*sin(0.5*pi*t);
        g=1+sum((x(2:V)-G).^2);
        f(2)=g*(1-(f(1)/g).^H);
    case 'DF2'
        G=abs(sin(0.5*pi*t));
        r=1+floor((n-1)*G);
        tmp=setdiff(1:n,r);
        g=1+sum((x(tmp)-G).^2);
        f(1)=x(r);
        f(2)=g*(1-(x(r)/g)^0.5);
    case 'DF3'
        G=sin(0.5*pi*t);
        H=G+1.5;
        g=1+sum((x(2:V)-G-x(1)^H).^2);
        f(1)=x(1);
        f(2)=g*(1-(x(1)/g)^H);
    case 'DF4'
        a=sin(0.5*pi*t);
        b=1+abs(cos(0.5*pi*t));
        c=max(abs(a), a+b);
        H=1.5+a;
        g=1+sum((x(2:V)-a*(x(1)/c).^2./[2:n]).^2);
        f(1)=g*abs(x(1)-a).^H;
        f(2)=g*abs(x(1)-a-b).^H;
    case 'DF5'
        G=sin(0.5*pi*t);
        w=floor(10*G);
        g=1+sum((x(2:end)-G).^2);
        f(1)=g*(x(1)+0.02*sin(w*pi*x(1)));
        f(2)=g*(1-x(1)+0.02*sin(w*pi*x(1)));
    case 'DF6'
       G=sin(0.5*pi*t);
        a=0.2+2.8*abs(G);
        y=x(2:V)-G;
        g=1+sum((abs(G)*y.^2-10*cos(2*pi*y)+10));
        f(1)=g*(x(1)+0.1*sin(3*pi*x(1))).^a;
        f(2)=g*(1-x(1)+0.1*sin(3*pi*x(1))).^a;
    case 'DF7'
        a=5*cos(0.5*pi*t);
        tmp=1/(1+exp(a*(x(1)-2.5)));
        g=1+sum((x(2:end)-tmp).^2);
        f(1)=g*(1+t)/x(1);
        f(2)=g*x(1)/(1+t);
        
    case 'DF8'
        G=sin(0.5*pi*t);
        a=2.25+2*cos(2*pi*t);
        tmp=G*sin(4*pi*x(1))/(1+abs(G));
        g=1+sum((x(2:end)-tmp).^2);
        f(1)=g*(x(1)+0.1*sin(3*pi*x(1)));
        f(2)=g*(1-x(1)+0.1*sin(3*pi*x(1))).^a;
        
    case 'DF9'
        N=1+floor(10*abs(sin(0.5*pi*t)));
        g=1;
        for i=2:n
            tmp=x(i)-cos(4*t+x(1)+x(i-1));
            g=g+tmp.^2;
        end
        f(1)=g*(x(1)+max(0, (0.1+0.5/N)*sin(2*N*pi*x(1))));
        f(2)=g*(1-x(1)+max(0, (0.1+0.5/N)*sin(2*N*pi*x(1))));
        
    case 'DF10'
        G=sin(0.5*pi*t);
        H=2.25+2*cos(0.5*pi*t);
        tmp=sin(4*pi*(x(1)+x(2)))/(1+abs(G));
        g=1+sum((x(3:end)-tmp).^2);
        f(1)=g*sin(0.5*pi*x(1)).^H;
        f(2)=g*sin(0.5*pi*x(2)).^H.*cos(0.5*pi*x(1)).^H;
        f(3)=g*cos(0.5*pi*x(2)).^H.*cos(0.5*pi*x(1)).^H;
    case 'DF11'
        G=abs(sin(0.5*pi*t));
        g=1+G+sum((x(3:end)-0.5*G*x(1)).^2);
        y=pi*G/6+(pi/2-pi*G/3)*x(1:2);
        f(1)=g*sin(y(1));
        f(2)=g*sin(y(2))*cos(y(1));
        f(3)=g*cos(y(2))*cos(y(1));
    case 'DF12'
        k=10*sin(pi*t);
        tmp1=x(3:end)-sin(t*x(1));
        tmp2=abs(sin(floor(k*(2*x(1:2)-1))*pi/2));
        g=1+sum(tmp1.^2)+prod(tmp2);
        f(1)=g*cos(0.5*pi*x(2)).*cos(0.5*pi*x(1));
        f(2)=g*sin(0.5*pi*x(2)).*cos(0.5*pi*x(1));
        f(3)=g*sin(0.5*pi*x(1));
    case 'DF13'
        G=sin(0.5*pi*t);
        p=floor(6*G);
        g=1+sum((x(3:end)-G).^2);
        f(1)=g*cos(0.5*pi*x(1)).^2;
        f(2)=g*cos(0.5*pi*x(2)).^2;
        f(3)=g*sin(0.5*pi*x(1)).^2+sin(0.5*pi*x(1)).*cos(p*pi*x(1)).^2 +...
            sin(0.5*pi*x(2)).^2+sin(0.5*pi*x(2)).*cos(p*pi*x(2)).^2;
    case 'DF14'
        G=sin(0.5*pi*t);
        g=1+sum((x(3:end)-G).^2);
        y=0.5+G*(x(1)-0.5);
        f(1)=g*(1-y+0.05*sin(6*pi*y));
        f(2)=g*(1-x(2)+0.05*sin(6*pi*x(2)))*(y+0.05*sin(6*pi*y));
        f(3)=g*(x(2)+0.05*sin(6*pi*x(2)))*(y+0.05*sin(6*pi*y));
    otherwise
        disp('no such test problem.')
end
end

