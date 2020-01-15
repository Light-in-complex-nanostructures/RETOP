function x=retexpm1(x,n);
% d'aprés Hervé Sauer 
% trés préférable à expm1 (plus rapide et au moins aussi précis)
% 
% Extension à d'autres degrés ( erreur relative < 1.e-11)
% y=retexpm1(x,n)
% si n=0   y=exp(x)
% si n=1   y=exp(x)-1
% si n=2   y=exp(x)-1-x
% si n=3   y=exp(x)-1-x-x^2/2
% si n=4   y=exp(x)-1-x-x^2/2-x^3/6
%
% See also EXPM1  LOG1P RETEXPM1  RETLOG1P


if nargin<2;n=1;end



switch n;
case 0;x=exp(x);% exp(x)
case 1;
[f,ff]=retfind(abs(x)>1);
x(f)=exp(x(f))-1;	
x(ff)=x(ff)/2;x(ff)=2*sinh(x(ff)).*exp(x(ff));% exp(x)-1
case 2;% exp(x)-1-x
[f,ff]=retfind(abs(x)>0.015);
x(f)=retexpm1(x(f),1)-x(f);
x(ff)=x(ff).^2.*polyval([1/720,1/120,1/24,1/6,.5],x(ff));
case 3;% exp(x)-1-x-x^2/3
[f,ff]=retfind(abs(x)>0.02);
x(f)=retexpm1(x(f),2)-.5*(x(f).^2);
x(ff)=x(ff).^3.*polyval([1/5040,1/720,1/120,1/24,1/6],x(ff));
case 4;% exp(x)-1-x-x^2/3-x^3/6
[f,ff]=retfind(abs(x)>0.05);
x(f)=retexpm1(x(f),3)-(x(f).^3)/6;
x(ff)=x(ff).^4.*polyval([1/40320,1/5040,1/720,1/120,1/24],x(ff));
end;	

% return;
% a=abs(x)>4.e-4;f=find(a);ff=find(~a);
% %[f,ff]=retfind(abs(x)>4.e-4);
% x(f)=exp(x(f))-1;
% x(ff)=x(ff).*polyval([1/24,1/6,.5,1],x(ff));


% retplot;for u=logspace(1,-20,100);x=(randn+i*randn)*u;xx=retprecis(x);y=exp(xx)-1;retplot(log10(u),log10([double(abs(y-expm1(x))),double(abs(y-retexpm1(x)))]));end;

% u=logspace(1,-20,100);x=(randn(size(u))+i*randn(size(u))).*u;xx=retprecis(x);y=exp(xx)-1;
% figure;font=retfont;loglog(u,double(abs(y-expm1(x))),'--r',u,double(abs(y-retexpm1(x))),'-k','linewidth',3);grid;xlabel('u=logspace(1,-20,100);x=(randn(size(u))+i*randn(size(u))).*u',font{:});set(legend('expm1','retexpm1'),font{:});set(gca,font{3:end})
% u=logspace(1,-20,100);x=randn(size(u)).*u;xx=retprecis(x);y=exp(xx)-1;
% figure;font=retfont;loglog(u,double(abs(y-expm1(x))),'--r',u,double(abs(y-retexpm1(x))),'-k','linewidth',3);grid;xlabel('u=logspace(1,-20,100);x=randn(size(u))*u',font{:});set(legend('expm1','retexpm1'),font{:});set(gca,font{3:end})
% n=1800;x=(randn(n)+i*randn(n))*1.e-2;tic;expm1(x);cpu_expm1=toc;tic;retexpm1(x);cpu_retexpm1=toc;cpu_expm1/cpu_retexpm1
% n=1800;x=randn(n)*1.e-2;tic;expm1(x);cpu_expm1=toc;tic;retexpm1(x);cpu_retexpm1=toc;cpu_expm1/cpu_retexpm1