function [mode,e,o,y,wy]=retpmode(pol,n_ext,n_int,L,k0,tracer,ordre);
%  [mode,e,o,y,wy]=retpmode(pol,n_ext,n_int,L,k0,tracer,ordre);
% 
% modes symétriques et antisymétriques  d'une couche de largeur L, d'indice n_int dans du milieu d'indice n_ext
%    Ce programme est spécialement adapté au modes d'une fente dans du metal 
%    mais fonctionne aussi pour d'autres configurations
%
% k0=2*pi/ld (1 par defaut)
% tracer 1 si on veut un tracé du champ du mode (0 par defaut)
% on peut aussi donner les parametres du tracer (voir rettchamp) 
% par exemple tracer=[1,13,i] tracera la partie reelle de E 
%
% la largeur suplementaire du tracé peut etre donnee avec L (L=[L,L_plus])
%
% si pol=0 E// ordre= 1 ,2 ,...
% ordre=  1 mode symetrique fondamental (par defaut) , 2m+1 ordres symetriques suivants
%         2 mode antisymetrique fondamental , 2m+2 autres ordres antisymetriques
% si pol=2 H// ordre= 0,1,  ...
% ordre=  0 mode symetrique fondamental (par defaut) , 2m ordres symetriques suivants
%         1 mode antisymetrique fondamental , 2m+1 autres ordres antisymetriques
% retourne une structure 'mode' de champs:
% 'constante_de_propagation'
% 'poynting' flux du vecteur du poynting,
% 'int_lorentz' integrale de lorentz 
%      du mode dont l'amplitude au centre de symetrie est donnée par le tableau suivant:
%                              
%                                              A  y
% ordre      sym     pol                          |
% 2*m+1        1      0       Ez=1   *******      |       *************
% 2*m+2       -1      0       Hy=1   *******   n_int      *** n_ext ***
% 2*m          1      2       Hz=1   *******      |       *************
% 2*m+1       -1      2       Ey=1   -------------+----------------------> x
% 	                                    <---- L ----->	
%
% on peut aussi sortir: e,o,y,wy
% le champ est alors tel que le flux du 'vecteur de Poynting complexe' .5*somme(E*Hx) soit égal à 1 (Lorentz=4) 
%
%%% EXEMPLES
%
% ld=.8;L=.1*ld; mode_S=retpmode(2,retindice(ld,2),1,L,2*pi/ld,1,0)
% ld=.8;L=.6*ld; mode_A=retpmode(2,retindice(ld,2),1,L,2*pi/ld,1,1)
%
% ld=.8;L=.1*ld; mode_S=retpmode(0,retindice(ld,2),1,L,2*pi/ld,1,1)
% ld=.8;L=.6*ld; mode_A=retpmode(0,retindice(ld,2),1,L,2*pi/ld,1,2)
%
% ld=.8;L=.06*ld; mode_A=retpmode(2,1,retindice(ld,2),[L,ld],2*pi/ld,[1,13,i],1)
%
% See also:RETPLASMON,RETMODE,RETTCHAMP
%



if nargin<7;ordre=0;end;
if nargin<6;tracer=[];end;
if ~isempty(tracer);if length(tracer)<=1;if tracer==0;tracer=[];else;tracer=[1:3,i];end;end;end;
if nargin<5;k0=1;end;
mu_int=n_int^pol;mu_ext=n_ext^pol;
sigma=pol-1;
sym=(pol==0& mod(ordre,2)==1) |(pol==2& mod(ordre,2)==0);

if length(L)==1;L_plus=L/2;else;L_plus=L(2);L=L(1);end;
% valeur de depart
particulier=0;
if imag(n_int)<.1*real(n_int) & imag(n_ext)<.1*real(n_ext) ;LL=inf; % tout dielectrique
else;% metal
plasmon=retplasmon(n_ext,n_int,2*pi/k0);
if imag(n_int)>.5*real(n_int);LL=0;else;LL=plasmon.penetration_dielectrique;end;% cas du metal à l'interieur
end;	

%if pol==2 & ((ordre==1 & L>1*plasmon.penetration_dielectrique)| (ordre==0 & L>2*plasmon.penetration_dielectrique)); %& plasmon.test_d_existence==1;% ¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤ modes plasmoniques
test_cad=0;
if pol==2 & ordre<2 & L>LL; %& plasmon.test_d_existence==1;% ¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤ modes plasmoniques
particulier=1;
gama=plasmon.constante_de_propagation;
for kk=1:50;
khi_int=khi(n_int,gama);
gama0=gama;

expifi=exp(i*khi_int*k0*L);
gama=retsqrt((2*khi_int^2*expifi./(expifi+(-1)^ordre*(1+expifi.^2)/2)*mu_ext^2-(mu_ext*n_int)^2+(mu_int*n_ext)^2)/(mu_int^2-mu_ext^2),-1);
if abs(gama0-gama)<.001*abs(gama) ;break;end;% test de convergence
end;
[gama,niter,er,erfonc,test_cad]=retcadilhac(@calf,struct('nitermin',5,'niter',50,'tol',1.e-10,'tolf',1.e-6),gama+[-i,0,i]*1.e-6,mu_int,mu_ext,n_int,n_ext,k0,L,sigma,ordre,particulier);

if all(test_cad);test=calf_test(gama,mu_int,mu_ext,n_int,n_ext,k0,L,sigma,ordre,particulier);test=abs(test);end;
end;
if ~all(test_cad);                                                                                     % ¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤ modes non plasmoniques

    fac=1;
	for kk=1:50; % iteration
	if sigma==1;
	A=fac*mu_int*2/(L*k0*mu_ext);	
	khi_int=roots([1,-2*ordre*pi/(k0*L),A^2+(ordre*pi/(k0*L))^2,0,(n_ext^2-n_int^2)*A^2]).';
	else;
	A=fac*mu_ext*2/(L*k0*mu_int);
	khi_int=roots([1,-2*ordre*pi/(k0*L),n_ext^2-n_int^2+(ordre*pi/(k0*L))^2+A^2,-2*ordre*pi/(k0*L)*(n_ext^2-n_int^2),(ordre*pi/(k0*L))^2*(n_ext^2-n_int^2)]).';
	end;
	khi_int=retsqrt(khi_int.^2,-1);
	khi_ext=retsqrt(khi_int.^2-n_int^2+n_ext^2,-1);
	st_warning=warning;warning off;prv=abs(khi_int*k0*L/2-ordre*pi/2+i*fac*(khi_ext.*mu_int./(khi_int.*mu_ext)).^sigma);warning(st_warning);
	f=find(prv<1.e-10);
	if length(f)<2;[prv,f]=sort(prv);f=f(1);end;
	khi_int=khi_int(f);khi_ext=khi_ext(f);
	talfa=-i*(khi_ext*mu_int./(khi_int*mu_ext)).^sigma;
	alfa=atan(talfa); % alfa./talfa,fac
	
	
	gama=retsqrt(n_int^2-khi_int.^2,-1);
    if length(gama)>1;
	test=calf(gama,mu_int,mu_ext,n_int,n_ext,k0,L,sigma,ordre,particulier);test=abs(test);
	[prv,ii]=min(test);% 'meilleure' solution
	gama=gama(ii);
	alfa=alfa(ii);
	talfa=talfa(ii);
    end;

	fac0=fac;%      fac  ,   gama

	if abs(talfa)<10*eps;fac=1;else;fac=.5*(fac+alfa./talfa);end;
	if isnan(fac);fac=fac0;end;
	%alfa
	if abs(fac0-fac)<.001*abs(fac) ;break;end;% test de convergence sur fac
	end;% kk

[gama,niter,er,erfonc,test]=retcadilhac(@calf,struct('nitermin',5,'niter',50,'tol',1.e-10),gama,mu_int,mu_ext,n_int,n_ext,k0,L,sigma,ordre,particulier);
test=calf(gama,mu_int,mu_ext,n_int,n_ext,k0,L,sigma,ordre,particulier);test=abs(test);
end;

if test>1.e-6;[gama,poynting,lorentz,khi_int,khi_ext,e,o,y,wy]=deal(nan+i*nan);%**************** mode non obtenu
else;                                                             %**************** mode obtenu
khi_int=khi(n_int,gama);khi_ext=khi(n_ext,gama);
if sym; % modes symetriques
poynting=real((gama*.25*L/mu_int)*(retsinc(i*L*imag(k0*khi_int)/pi)+retsinc(L*real(k0*khi_int)/pi)))...
	                         +.5*abs(cos(k0*khi_int*L/2))^2*real(gama/mu_ext)/imag(khi_ext*k0);
lorentz=gama*L/mu_int*(1+retsinc(k0*L*khi_int/pi))+2i*gama*cos(k0*khi_int*L/2)^2/(mu_ext*khi_ext*k0);
else;   % modes antisymetriques
poynting=real((gama*.25*L/mu_int)*real(retsinc(i*L*imag(k0*khi_int)/pi)-retsinc(L*real(k0*khi_int)/pi)))...
	                         +.5*abs(sin(k0*khi_int*L/2))^2*real(gama/mu_ext)/imag(khi_ext*k0);
lorentz=gama*L/mu_int*(1-retsinc(k0*L*khi_int/pi))+2i*gama*sin(k0*khi_int*L/2)^2/(mu_ext*khi_ext*k0);
% pour avoir Hy=1 en y=0
poynting=poynting*abs(mu_int/khi_int)^2;
lorentz=lorentz*(mu_int/khi_int)^2;
end;  % modes symetriques ?
end;                                                              %**************** mode obtenu ?

mode=struct('constante_de_propagation',gama,'poynting',poynting,'int_lorentz',lorentz,'khi_int',khi_int,'khi_ext',khi_ext,'test',test);

if (~isempty(tracer) | nargout>1) & ~isnan(gama);
[y,wy]=retgauss(0,L/2+L_plus,20,1+round(5*(L/2+L_plus)*k0/(2*pi)),L/2);[f_int,f_ext]=retfind(y<L/2);	
e=zeros(length(y),1,4);o=zeros(length(y),1,1);
o(f_int,:,:)=k0*n_int;o(f_ext,:,:)=k0*n_ext;o=[flipdim(o,1);o];
if sym; % modes symetriques
e(f_int,:,1)=cos(khi_int*k0*y(f_int));
e(f_int,:,2)=-khi_int/mu_int*sin(khi_int*k0*y(f_int));
e(f_int,:,3)=-i*gama/mu_int*cos(khi_int*k0*y(f_int));
e(f_int,:,4)=-k0*khi_int*sin(khi_int*k0*y(f_int));

e(f_ext,:,1)=cos(khi_int*k0*L/2)*exp(i*khi_ext*k0*(y(f_ext)-L/2));
e(f_ext,:,2)=i*khi_ext/mu_ext*cos(khi_int*k0*L/2)*exp(i*khi_ext*k0*(y(f_ext)-L/2));
e(f_ext,:,3)=-i*gama/mu_ext*cos(khi_int*k0*L/2)*exp(i*khi_ext*k0*(y(f_ext)-L/2));
e(f_ext,:,4)=i*k0*khi_ext*cos(khi_int*k0*L/2)*exp(i*khi_ext*k0*(y(f_ext)-L/2));

y=[-fliplr(y),y];wy=[fliplr(wy),wy];
e=[flipdim(e,1);e];
f=find(y<0);
e(f,:,[2,4])=-e(f,:,[2,4]);
else;                                                     % modes antisymetriques
e(f_int,:,1)=mu_int/khi_int*sin(khi_int*k0*y(f_int));
e(f_int,:,2)=cos(khi_int*k0*y(f_int));
e(f_int,:,3)=-i*gama/khi_int*sin(khi_int*k0*y(f_int));
e(f_int,:,4)=i*k0*gama*sin(khi_int*k0*y(f_int));

e(f_ext,:,1)=mu_int/khi_int*sin(khi_int*k0*L/2)*exp(i*khi_ext*k0*(y(f_ext)-L/2));
e(f_ext,:,2)=i*khi_ext/mu_ext*mu_int/khi_int*sin(khi_int*k0*L/2)*exp(i*khi_ext*k0*(y(f_ext)-L/2));
e(f_ext,:,3)=-i*gama/mu_ext*mu_int/khi_int*sin(khi_int*k0*L/2)*exp(i*khi_ext*k0*(y(f_ext)-L/2));
e(f_ext,:,4)=i*k0*khi_ext*mu_int/khi_int*sin(khi_int*k0*L/2)*exp(i*khi_ext*k0*(y(f_ext)-L/2));

y=[-fliplr(y),y];wy=[fliplr(wy),wy];
e=[flipdim(e,1);e];
f=find(y<0);
e(f,:,[1,3,4])=-e(f,:,[1,3,4]);

e=i*e; % pour avoir,aprés declonage Hy=1 en y=0 et pas -1

end;	
e(:,:,2:4)=-i*e(:,:,2:4); % declonage
e=2*e/sqrt(lorentz);
end;
if ~isempty(tracer)& ~isnan(gama);x=0;rettchamp(e,o,x,y,pol,tracer,[],[],rettexte(gama,niter,n_int,n_ext,L,k0));end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function f=calf_test(gama,mu_int,mu_ext,n_int,n_ext,k0,L,sigma,ordre,particulier);% modif 2010
khi_int=khi(n_int,gama);khi_ext=khi(n_ext,gama);
f=tan(k0*L*khi_int/2+ordre*pi/2)+i*((khi_ext./khi_int)*(mu_int/mu_ext)).^sigma;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function f=calf1(gama,mu_int,mu_ext,n_int,n_ext,k0,L,sigma,ordre,particulier);
khi_int=khi(n_int,gama);khi_ext=khi(n_ext,gama);
alfa=atan(-i*((khi_ext./khi_int)*(mu_int/mu_ext)).^sigma);
f=(k0*L)*khi_int-(ordre*pi+2*alfa);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function f=calf(gama,mu_int,mu_ext,n_int,n_ext,k0,L,sigma,ordre,particulier)
if ordre>1 | ~particulier;f=calf1(gama,mu_int,mu_ext,n_int,n_ext,k0,L,sigma,ordre);return;end;% seul le premier ordre peut 'zapper'
khi_int=khi(n_int,gama);
expifi=exp(i*khi_int*k0*L);
if sigma==1;% pol=2
f=2*(n_int.^2-gama.^2).*expifi./(expifi+(-1)^ordre*(1+expifi.^2)/2)-(gama.^2*(mu_int^2-mu_ext^2)+(mu_ext*n_int)^2-(mu_int*n_ext)^2)/(mu_ext^2);
else; % pol=0
f=2*(n_int.^2-gama.^2).*expifi./(expifi-(-1)^ordre*(1+expifi.^2)/2)-n_int^2+n_ext^2;
end;	

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function kkhi=khi(n,gama); 
kkhi=retsqrt(n^2-gama.^2,-1);


