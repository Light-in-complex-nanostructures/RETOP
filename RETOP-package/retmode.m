function varargout=retmode(varargin)
%%%%%%%%%%%%%%%%%%%%%%%
%%%   version 1 D    %%
%%%%%%%%%%%%%%%%%%%%%%%
%
%  [gama,kmax,vm,er,e,o,y,wy,vmc]=retmode(pol,n,h,ga,gb,hb,zs,k0);
%
%  milieux stratifiés (0D) recherche des modes(ou tracé préliminaire pour detecter la présence de modes);
%  pol:polarisation (0  E   2 H)
%
%  Si milieux homogènes en haut et en bas:
%    n:tableau des indices  h:tableau des epaisseurs ( de haut(z grand) en bas( z petit)  length(n)=length(h)+2)
%   (on peut aussi coder des couches de graphéne: h=i,n=sig_G) 8 2012
%
%  Si modes de Bloch en haut ou en bas:
%    n={nh,nc,nb}: nh,nc,nb tableau des indices en haut, au centre, en bas 
%    h={hh,hc,hb}: hh,hc,hb tableau des epaisseurs en haut, au centre, en bas 
%    length(nh)=length(hh),length(nc)=length(hc),length(nb)=length(hb)
%    Si d'un coté on a un milieu homogene,on peut prendre hh (ou hb) nul 0 ou vide 
%	Remarque : les conditions aux limites pour les modes de Bloch (attenuation)  ne sont pas les mêmes que pour les milieux homogènes (coupure de Maystre)
%  Dans le cas de milieux homogènes, le programme le reconnaît et revient à la convention correspondante
%
%  ga: tableau des valeurs initiales des modes 
% (si absent tracé de r et t de gb(1) a gb(2) en gb(3) points  (1000 par defaut) alors gama est le tableau des max qui peut servir de point de depart pour la recherche precise des modes)
%  gb:bornes et nombre de points pour le tracé fonction de beta (facultatif)
%  hb([1,2]):hauteurs en haut et en bas pour le tracé du champ des modes (facultatif: si absent pas de tracé)
%     ( si modes de Bloch en haut ou en bas le tracé est arrondi à un nombre entier de périodes (au moins 1) )
%  si hb a 3 termes:tracé de l'intensité, si hb a seulement 2 termes tracé de real et imag du champ
%  gama: modes trouves  kmax:nombre de pas (meme dimension que ga) si kmax=100 le mode n'a pas été trouvé
%
%  vm,er volume modal et champ sur la source placée en xs, normalisé
%-------------------------------------------------------
% vm(1)=somme(abs(eps*E^2)) /max (abs(eps* abs(E)^2);
% vm(2)=somme(abs( mu*H^2)) /max  (abs(mu* abs(H)^2);
% vm(3)=somme(abs(eps*E^2)+abs(mu*H^2))/max( abs(eps* E^2)+abs(mu* H^2));
% er:  composantes du champ sur la source (zs) normalises par le volume reel=somme(abs(eps*E^2)+abs(mu*H^2))
%
%  zs position de la source depuis le bas (ne sert que pour le calcul de er)
%   ( Si modes de Bloch en haut ou en bas zs est corrigé de l'arrondi du au nombre entier de periodes en bas)
%  si on ne veut pas le tracé (mais le calcul des champs) hb=[hb(1),hb(2),0]
%
%  le calcul du vm et er est PROLONGE A L'INFINI en haut et en bas en utilisant la décroissance exponentielle des champs
%   (le resultat est donc indépendent de hb)
%  Les max de eps*abs(E)^2,mu*abs(H)^2,eps*abs(E)^2+mu*abs(H)^2 correspondent au domaine où les champs sont calculés  
%
% e,o,y,wy: si calcul des champs. Ce sont des cell-array dont la dimension est celle de gama, ce qui permet de tracer les champs avec rettchamp.
% 
% vmc=somme(epsz*E^2+mux*Hx^2-muy*Hy^2)/(.5*somme(E*Hx))   volume modal complexe 
%   (la vitesse de groupe 'géométrique' est vg=4*k0/vmc)
% le champ est tel que le flux du 'vecteur de Poynting complexe' .5*somme(E*Hx) soit égal à 1 (Lorentz=4) 
%
% Balayage dans le plan complexe, recommandé pour une recherche systematique
%.....................................................................................
% [zz,kmax,z]=retmode(pol,n,h,betr,beti,k0,hb);
% pol,n,h,: comme retmode 1D
% betr,beti vecteurs des neff complexes (K//=(Betr+i*Beti)*k0, où [Beti,Betr]=ndgrid(beti,betr);
% betr et beti doivent avoit au moins 2 elements 
% hb si on veut le tracé des modes
% Le programme calcule t avec retabeles vectorialisé et affiche une quantité deduite de t par utilisation du theoreme de Cauchy, dont les maximas correspondent aux poles.
% On peut ensuite selectionner des points de depart z
% si zz est en sortie,les modes sont cherchés (valeur zz et kmax nombre d'iterations)
% .....................................................................................
% Suivi de mode
% [lld,ggama,coupure]=retmode(pol,indice,h,hb,ld1,gama1,parm);
% indice:fonction	 n=indice(ld) est le parametre de retmode
% pol,h, parametres de retmode
% hb
% ld1,gama1=K//1=neff0*k0 point de départ approximatif
% le mode de depart est cherché et eventuellement tracé (selon hb)
% 
% parm=struct('gama_reel',0,'antizap',1.5,'tracer',0, ...
% 	'ldmin',-inf,'ldmax',inf,'dld0',1.e-6,'dldmax',5.e-3,'dldmin',1.e-11, ...
% 	'gamamin',-inf,'gamamax',inf,'dgama0',1.e-6,'dgamamax',5.e-3,'dgamamin',1.e-11);
%                       par defaut
% gama_reel=0 on pole en gama=K a ld reel fixé Indice est une fonction complexe de variable reelle
% gama_reel=1 on pole en ld a gama reel fixé Indice est une fonction complexe de variable complexe (modele de type Drude necessaire)
% ldmin, ldmax bornes en ld ( ou gamamin, gamamax bornes en gama)
% dld0 (ou dgama0)pas initial (ou dgama0)
% dldmax (ou dgamamax)  plus grand pas permis
% dldmin (ou dgamamin)  plus petit pas permis ( coupure)
% Le mode est suivi dans les 2 sens depuis ld1 (ou gama1)
% au cours du calcul, le pas peut être augmenté (si kmax<5) jusqu'à dldmax (ou dgamamax)
% Quand on s'approche des coupures de mode, le pas diminue jusqu'à dldmin (ou dgamamin)
% Si la variation est trop brusque ( parametre antizap>1) , on diminue le pas
%  en sortie:
% lld longueurs d'ondes  ggama= K//=neff*k0
% coupure vecteur de 2 elements: si le mode arrive à un cut_off à gauche coupure(1)=1, si le mode arrive à un cut_off à droite coupure(2)=1
% si indice donne un cellarray, c'est en fait epsilon et on utilise retyeh
% 
% %%  exemples 1D
%%%**************
%
% retmode(0,[1,3.5,1],1) ;                             % trace de rh rb et t fonction de gama de 1 a 3.5
% retmode(0,[1,3.5,1],1,[],[-1,1,101]) ;               % trace de rh rb et t fonction de gama de -1 a 1 avec 101 points
% [gama,kmax]=retmode(0,[1,3.5,1],1,[3,1.1])           % recherche du mode voisin de 3 
% [gama,kmax]=retmode(0,[1,3.5,1],1,[3,1.1],[],[2,2])  % recherche et trace des modes voisins de 3 et 1.1 (real imag E )
% [gama,kmax]=retmode(0,[1,3.5,1],1,3,[],[2,2,1])      % recherche et trace du mode voisin de 3  (abs(E)^2)
% [gama,kmax,vm,erm]=retmode(0,[1,3.5,1],1,[3,1.1],[],[2,2],2.5) 
%%% aussi calcul du volume modal et du champ normalise(source en 2.5)
%
% ld=.97;k0=2*pi/ld;[gama,kmax]=retmode(2,[1,retindice(ld,1)],[],[1.01*k0],[],[1,1,.1],1,k0)          % plasmon  %
%
%  n2=2;n3=3;h2=2*pi/(4*n2);h3=2*pi/(4*n3);h0=4*h3;
%  mh=6;mb=8;h=[repmat([h2,h3],1,mh-1),h2,h0,repmat([h2,h3],1,mb-1),h2];n=[n3,repmat([n2,n3],1,mh+mb)];
% [gama,kmax,vm,erm]=retmode(0,n,h,[0.02i],[],[h3,h3],mb*(h2+h3)+h0/2)
%
%%% Exemple avec modes de Bloch en haut ou en bas
% [gama,kmax]=retmode(2,{[3,.1+5i],[],1},{[.1,.1],[],0},[1.2],[1,2.5],[1,2,1]);
% ld=.8;k0=2*pi/ld;[gama,kmax]=retmode(0,{[1,1.5],3,[2,1.4]},{[.1,.2],.2,[.2,.3]},[1.4,1.9,2.7]*k0,[],[1,1,1],1.1,k0);
% 
%  n2=2;n3=3;h2=2*pi/(4*n2);h3=2*pi/(4*n3);h0=4*h3;[gama,kmax,vm,er,e,o,y,wy,vmc]=retmode(0,{[n3,n2],n3,[n2,n3]},{[h3,h2],h0*1.00001,[h2,h3]},[.008],[],30*[h3,h3,1]);
%  retpoynting(e{1},[1,0],[],wy{1})
%
%%% Balayage
% ld=.54;k0=2*pi/ld;nm=retindice(ld,1.7);[neff,kmax]=retmode(2,[1,nm,3.5,nm],[.025,.1],linspace(-8,16,400),linspace(-1,17,400),k0,[.1,.1,1]);
%%% ld=.54;k0=2*pi/ld;nm=retindice(ld,1.7);for teta=linspace(.6*pi,.1*pi,50);retsqrt(0,0,teta);retmode(2,[1,nm,3.5,nm],[.025,.1],linspace(-8,16,400),linspace(-1,17,400),k0,[.1,.1,1]);title(rettexte(teta));drawnow;close gcf;end;retsqrt(0,0,pi/2);
%%% La coupure varie
%%% Suivi
% [lld,ggama]=retmode(2,@(ld) [1,1.5,1.2],.5,[.5,.5,1],.6,1.22*(2*pi/.6),struct('ldmin',.5,'ldmax',1));figure;neff=ggama.*lld/(2*pi);plot(lld,real(neff),'.-k')
%%% recherche automatique des modes TM ( attention à bien fixer les limites du balayage initial)
% k0=2*pi/.6;nh=1;n1=1.5;nb=1.2;nm=0.1+5i;
% gamaH=retmode(2,[nh,n1,nm,nb],[.3,.05],[],k0*[.9*nh,1.5*n1],[],[],k0);close gcf;
% [gamaH,kmax]=retmode(2,[nh,n1,nm,nb],[.3,.05],gamaH,[],[],[],k0);gamaH=retelimine(gamaH(kmax<100 & abs(imag(gamaH)<.1*k0)));
% [gamaH,kmax]=retmode(2,[nh,n1,nm,nb],[.3,.05],gamaH,[],[.1,.1],[],k0);% pour le tracé des modes selectionnes
%
%    %%%%%%%%%%%%%%%%%%%%%
%    %    version 2 D    %
%    %%%%%%%%%%%%%%%%%%%%%
%
% [gama,kmax,vm,erm,ept,e,o,x,z,wx,wz,vmc]=retmode(init,mm,pols,u,ga,gb,hb,fich);
% modes d'un objet cylindrique de section donnée par recherche de pôle en beta conique
% l'objet est décrit par un 'maillage de texture'  u=retu(d,... en 2 D
%
%
%
%  ^
%  |
%  |              ****
%  |           *************
%  |  d(2)    ****************
%  |        *******   ********
%  |        ******  S  ********-   < - - -source au centre de symetrie:sym(3:4);
%  |       *******   ********                (même si pas de symetrie) 
%  |           ************          
%  |             ********       
%  |                **            
%  |                            X
%  |------------------------------->  Y=cao(2) (même si pas de cao)
%                   d(1)
%
%
% mm:ordres de fourier (de -mm a mm) ou mm(1) a mm(2)
% pour des raisons de stabilite numérique (surtout pour le calcul des champs), 
% il est conseillé de travailler en G ( mm=5+i)
% init crée par retinit 2 D 
%
%  ga tableau des valeurs initiales des modes
%   (si absent ou vide trace de la reponse de la source pour les valeurs de pols en 50 points en beta entre gb(1) et gb(2) )
%  gb:vecteur de longueur 3 :bornes pour le tracé fonction de beta et nombre de points (facultatif) (gb(3)=20 par defaut)
%  hb:vecteur de longueur 2  hauteurs en bas et en haut pour le tracé du champ des modes  (facultatif: si absent pas de tracé)
%  si on ne veut pas le tracé (mais le calcul des champs) hb=[hb(1),hb(2),0]
%        fich nom de fichier (chaine de caracteres) pour stocker les resultats intermediaires:gama z zz er k (facultatif)
%        cette utilisation de fich est obsolete,on peut aussi au lieu de fich donner les parametres de retcadilhac
%        struct('niter',30,'nitermin',0,'tol',eps,'tolf',inf) par defaut: struct('niter',30,'tol',1.e-10,'tolf',1.e-5)
%  au retour :gama modes trouves  kmax:nombre de pas (meme dimension que ga)
%
% vm,erm volume modal et champ sur la source normalisé
%-------------------------------------------------------
% vm(1)=somme(eps*abs(E)^2)/(max eps* E^2);
% vm(2)=somme(mu*abs(H)^2)/(max  mu* H^2);
% vm(3)=somme(eps*abs(E)^2+mu*abs(H)^2)/max( eps* E^2+mu* H^2);
% er:  composantes du champ sur la source normalise par le volume reel
% vm(1)=somme(abs(eps*E^2)) /max (abs(eps* abs(E)^2);
% vm(2)=somme(abs( mu*H^2)) /max  (abs(mu* abs(H)^2);
% vm(3)=somme(abs(eps*E^2)+abs(mu*H^2))/max( abs(eps* E^2)+abs(mu* H^2));%    erm:  composantes du champ sur la source normalise par le volume reel
% ept:  composantes du champ sur la source normalisées par le flux du  vecteur de poynting
%  ces composantes sont dans la base initiale 2D où l'objet est décrit par le  'maillage de texture'
%  si hb n'est pas indiqué on prend la hauteur de l'objet depuis la source de chaque cote ( pas de trace)
%  si on ne veut pas le trace hb=[hb(1),hb(2),0]
%  sur les tracés de champs,la source est indiquée par un rectangle blanc
% 
% vmc=somme(epsx*Ex^2-epsy*Ey^2+epsz*Ez^2+mux*Hx^2-muy*Hy^2+muz*Hz^2)/(.5*somme(Ez*Hx-Ex*Hz))   volume modal complexe 
%   (la vitesse de groupe 'géométrique' est vg=4*k0/vmc)
% le champ est tel que le flux du 'vecteur de Poynting complexe' .5*somme(Ez*Hx-Ex*Hz) soit égal à 1  
%
%
%%  EXEMPLES 2D
%%**************
% ld=2*pi;k0=2*pi/ld;d=[8,5];m=5;lcao=[2,2];sym=[1,1,0,0];cao=[-d/2,lcao];init=retinit(d,[0,0,0,0],[],sym,cao);
% u=retu(init,{1,[0,0,5,4,3,5],[0,0,1.5,1,1,1],k0});
% pols=[1:6];retmode(init,m+i,pols,u,[],[1.1,3]); % trace de la reponse de la source 
% pols=6;[gama,kmax]=retmode(init,m+i,pols,u,[2.5],[],[1,1]),% recherche des modes et trace
%
% [uu,dd]=retrotation(1,u,d);ccao=[dd/2,lcao];ssym=sym;ssym([1,2])=-sym([2,1]);
% iinit=retinit(dd,[0,0,0,0],[],ssym,ccao);
% pols=3;[gama,kmax]=retmode(iinit,m+i,pols,uu,[1.5],[],[1,1])
%
% d=[8,5];m=5;lcao=[2,2];sym=[-1,1,0,0];cao=[d/2,lcao];init=retinit(d,[0,0,0,0],[],sym,cao);
% u=retu(init,{1,[0,0,5,4,3,5],[0,0,1.5,1,1,1],k0});
% pols=[1:6];retmode(init,m+i,pols,u,[],[1.1,3]); % trace de la reponse de la source 
% pols=[2,3];[gama,kmax]=retmode(init,m+i,pols,u,[2],[],[1,1]),% recherche des modes et trace
%
% % See also RETPMODE RETVM RETCADILHAC RETMARCUSE RETINDICE RETYEH


%%%%%%%%%%%%%%%%%%%%
%   AIGUILLAGE     %
%%%%%%%%%%%%%%%%%%%%




if iscell(varargin{1});[varargout{1:nargout}]=ret2mode(nargout>2,varargin{:});return;
else;
if iscell(varargin{3});
[varargout{1:nargout}]=ret1mode_bloch(nargout>2,varargin{:});return;
else;
if ~isnumeric(varargin{2});[varargout{1:nargout}]=suit_mode(varargin{:});return;end;
if nargin>5 & length(varargin{6})==1;[varargout{1:nargout}]=balayage_1D(varargin{:});return;
else; [varargout{1:nargout}]=ret1mode(nargout>2,varargin{:});return;end;
end;
end;




% if iscell(varargin{1});[varargout{1:nargout}]=ret2mode(nargout>2,varargin{:});% retmode 2D
% elseif iscell(varargin{2});[varargout{1:nargout}]=ret1mode_bloch(nargout>2,varargin{:});% mode 1D avec modes de bloch
% elseif nargin>5 & length(varargin{6})==1 [varargout{1:nargout}]=balayage_1D(varargin{:});% balayage
% elseif isnumeric(varargin{2}) [varargout{1:nargout}]=ret1mode(nargout>2,varargin{:});% retmode 1D
% else [varargout{1:nargout}]=suit_mode(nargout>2,varargin{:});% suivi de mode
% end;

%%%%%%%%%%%%%%%%%%%%
%   version 1 D    %
%%%%%%%%%%%%%%%%%%%%


function [gama,kmax,vm,erm,e,o,yy,wyy,vmc]=ret1mode(calvm,pol,n,h,ga,gb,h_b,ys,k0);
vm=[];erm=[];ept=[];vmc=[];
if nargin<5;ga=[];end;if nargin<6;gb=[];end;if nargin<7;h_b=[];end;if nargin<8;ys=0;end;if nargin<9;k0=1;end; 
if isempty(ys);ys=0;end;
tracer=(nargin>=7)&(~isempty(h_b));if (length(h_b)==3)&(h_b(3)==0) tracer=0;end;
if isempty(h_b);h_b=[0,0];end;
[gama,vmc]=deal(zeros(size(ga)));kmax=gama;[vm,erm,ept]=deal(zeros(3,length(ga)));[e,o,yy,wyy]=deal(cell(1,length(ga)));

if isempty(ga); % trace de rh rb t
figure;
if length(gb)==3;nn=gb(3);else nn=1000;end;
if isempty(gb);x=linspace(min([real(n(1)*k0),real(n(end)*k0)]),max(real(n*k0)),nn);else;x=linspace(gb(1),gb(2),nn);end;

tab=[[0;h(:);0],n(:)];init={pol,x,k0};sh=retb(init,n(1),1); sb=retb(init,n(end),-1); 
s=reteval(retss(sh,retcouche(init,tab),sb));
rh=squeeze(s(1,2,1,:));rb=squeeze(s(2,1,1,:));t=squeeze(s(1,1,1,:));
%rh=zeros(1,nn);rb=zeros(1,nn);t=zeros(1,nn); for ii=1:nn;[rh(ii),rb(ii),prv,t(ii)]=calg(x(ii),n,k0,h,pol+1i);end;t=1./t;
plot(x,abs(rh).^2,'-k',x,abs(rb).^2,'--k',x,abs(t).^2,':k','linewidth',2);legend('rh','rb','t');xlabel(['g    deltag=',num2str(x(end)-x(1)),' g= ',num2str(.5*(x(1)+x(end)),10)]);
% recherche des max
gama=[];
if max(abs(rh))>10*eps;[xe,ye,max_ou_min]=retextremum(x,abs(rh).^2,nan);gama=[gama,xe(max_ou_min==1)];end;% pour eliminer le guide homogene ...
if max(abs(rh))>10*eps;[xe,ye,max_ou_min]=retextremum(x,abs(rb).^2,nan);gama=[gama,xe(max_ou_min==1)];end;
[xe,ye,max_ou_min]=retextremum(x,abs(t).^2,nan);gama=[gama,xe(max_ou_min==1)];gama=retelimine(gama);

else % recherche des modes 
[gama,kmax]=deal(zeros(size(ga)));
for ii=1:length(ga);
for kkk=0:3;% modif 2015
[gama(ii),kmax(ii),prv,prv,test]=retcadilhac(@calg,struct('nout',4,'niter',100,'nitermin',0,'tol',1.e-10,'tolf',1.e-5),ga(ii),n,k0,h,pol+1i*kkk);
if all(test);break;end;
end;   
end
end;

% trace du champ des modes et calcul du volume modal
%tracer=(nargin>=7)&(~isempty(h_b));if (length(h_b)==3)&(h_b(3)==0) tracer=0;end;
if tracer;figure;end;
if calvm|tracer;% *************************

kkk=ceil(sqrt(length(ga)));kkkk=0;
for ii=1:length(ga);%%%%%%%%%%%%%%
[rh,rb,t,gg,e{ii},yy{ii},wyy{ii},o{ii},khih,khib,vmc(ii)]=calg(gama(ii)*(1+1.e-14),n,k0,h,pol,h_b,ys);% *(1+1.e-14) modif 19 11 2012 bugg Anthony
 
if calvm;xx=0;wxx=1;
[prv,erm(:,ii),ec,ener,enere,enerh,poynting]=retvm(e{ii},o{ii},{xx,yy{ii}},{wxx,wyy{ii}},[0,ys]);vm(:,ii)=prv(:);
% termes correctifs dus au prolongement exponentiel decroissant à l'infini

exph=2*imag(khih);expb=2*imag(khib);
eneree=enere+(abs(e{ii}(end,1,1)).^2).*abs(o{ii}(end,1,1))/exph+(abs(e{ii}(1,1,1)).^2).*abs(o{ii}(1,1,1))/expb;
enerhh=enerh+(abs(e{ii}(end,1,2)).^2).*abs(o{ii}(end,1,2))/exph+(abs(e{ii}(1,1,2)).^2).*abs(o{ii}(1,1,2))/expb...
+(abs(e{ii}(end,1,3)).^2).*abs(o{ii}(end,1,3))/exph+(abs(e{ii}(1,1,3)).^2).*abs(o{ii}(1,1,3))/expb;

erm(:,ii)=erm(:,ii)*sqrt((enere+enerh)/(eneree+enerhh));
vm(1,ii)=vm(1,ii)*(eneree/enere);
vm(2,ii)=vm(2,ii)*(enerhh/enerh);
vm(3,ii)=vm(3,ii)*(eneree+enerhh)/(enere+enerh);

else vm=[];erm=[];end;

if tracer;% ..................
kkkk=kkkk+1;retsubplot(kkk,kkk,kkkk);
if length(h_b)==2;
plot(yy{ii},real(e{ii}(:,1,1)),'-k',yy{ii},imag(e{ii}(:,1,1)),'--k','linewidth',2);%if kkkk==1;set(legend('real','imag'),'fontsize',8);end;
else;
plot(yy{ii},abs(e{ii}(:,1,1)).^2,'-k','linewidth',2);%if kkkk==1;set(legend('intensite'),'fontsize',8);end;
end;
% %set(leg,'fontsize',max(7,7-kkk));
% set(leg,'fontsize',8);
%xlabel([' gama=',num2str(gama(ii),10),' kmax=',num2str(kmax(ii))],'fontsize',max(7,7-kkk));
%xlabel([' gama=',num2str(gama(ii),10),' kmax=',num2str(kmax(ii))],'fontsize',8);
xlabel(['neff=',num2str(gama(ii)/k0,10),' kmax=',num2str(kmax(ii))],'fontsize',10);set(gca,'YTick',[]);
hold on;axis tight;ax=axis;oo=real(sqrt((o{ii}(:,1,1).*o{ii}(:,1,2))));oo=ax(3)+.5*oo*(ax(4)-ax(3))/max(abs(oo));plot(yy{ii},real(oo),':k','linewidth',2);hold off;axis tight;
retfont(gcf,0);
end  % ..................

end;  %%%%%%%%%%%%%%

end;% *************************

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [rh,rb,t,gg,e,yy,wyy,o,khih,khib,vmc]=calg(g,n,k0,h,pol,h_b,ys);
if iscell(n);% yeh n, est eps
[rh,rb,t,e,yy,wyy,o,khih,khib,vmc]=deal(nan);
gg=calg_yeh(k0,n,g,h,pol);
return;
end;
	
parm=imag(pol);pol=real(pol);
n=n(:);h=h(:);
% graphéne
[f_pas_graphene,f_graphene]=retfind(imag(h)==0);
if ~isempty(f_graphene);
sig_G=n(f_graphene+1);
n=n([1;f_pas_graphene+1;end]);h=h(f_pas_graphene);
f_graphene=f_graphene-[0:(length(f_graphene)-1)].';% numero des dioptres à graphéne
end;


y=-[0;cumsum(h)];
ep=retep(n,pol,k0);mu=ep(2,:).';ep=ep(1,:).';
m=length(n);
khi=retbidouille(retsqrt(ep.*mu-g^2,-1));
y=[y(1);y;y(end)];
p=(1:m-1).';
prv1=exp(-i*khi(p).*(y(p+1)-y(p)));
prv2=exp(i*khi(p+1).*(y(p+1)-y(p+2)));
prv3=khi(p).*mu(p+1)./(mu(p).*khi(p+1));
[ii,jj,M]=deal([]);
ii=[ii;p];jj=[jj;p];M=[M;ones(m-1,1)];
ii=[ii;p];jj=[jj;p+m];M=[M;prv1];
ii=[ii;p];jj=[jj;p+m+1];M=[M;-ones(m-1,1)];
ii=[ii;p];jj=[jj;p+1];M=[M;-prv2];

ii=[ii;p+m-1];jj=[jj;p];M=[M;prv3];
ii=[ii;p+m-1];jj=[jj;p+m];M=[M;-prv1.*prv3];
ii=[ii;p+m-1];jj=[jj;p+m+1];M=[M;ones(m-1,1)];
ii=[ii;p+m-1];jj=[jj;p+1];M=[M;-prv2];

% graphéne
if ~isempty(f_graphene);
if pol==0;% TE
ii=[ii;f_graphene+m-1];jj=[jj;f_graphene+1];M=[M;-1i*sig_G.*prv2(f_graphene).*mu(f_graphene+1)./khi(f_graphene+1)];    
ii=[ii;f_graphene+m-1];jj=[jj;f_graphene+m+1];M=[M;-1i*sig_G.*mu(f_graphene+1)./khi(f_graphene+1)];    
else;% TM
ii=[ii;f_graphene];jj=[jj;f_graphene+1];M=[M;-1i*sig_G.*prv2(f_graphene).*khi(f_graphene+1)./mu(f_graphene+1)];    
ii=[ii;f_graphene];jj=[jj;f_graphene+m+1];M=[M;1i*sig_G.*khi(f_graphene+1)./mu(f_graphene+1)];    
end;
end;

M=sparse(ii,jj,M,2*m-2,2*m);
if nargout<5;% <****   calcul de r et t
a=-M(:,[1:m-1,m+2:2*m])\M(:,[m,m+1]);
rh=a(1,2);
rb=a(2*m-2,1);
t=a(1,1);
% gg=1/t;
switch parm;% modif 1 2015
case 0;gg=1/t;
case 1;
try;gg=eigs(M(:,[1:m-1,m+2:2*m]),1,'SM',struct('disp',0));catch; gg=nan;end;
case 2;gg=1/rh;
case 3;gg=1/rb;
end;

else;        %  <****   calcul du champ
	
[rh,rb,t,gg,yy,wyy]=deal([]);

st_warning=warning;warning off;try;[a,prv]=eigs(M(:,[1:m-1,m+2:2*m]),1,'sm',struct('disp',0));catch;a=-M(:,[1:m-1,m+2:2*m])\[1;zeros(2*m-3,1)];end;
a=reshape([a(1:m-1);0;0;a(m:2*m-2)],m,2);

y(1)=y(1)+h_b(1);y(m+1)=y(m+1)-h_b(2);y=y-y(end);
e=zeros(0,1,3);o=zeros(0,1,3);
for ii=m:-1:1;
	
[yyy,wyyy]=retgauss(y(ii+1),y(ii),20,-1*min(100,ceil(eps+4*abs(khi(ii))*abs(y(ii)-y(ii+1)))));% eps ajout 9 2012 pour eviter une fausse normalisation si hb(2)=0 ou hb(1)=0
if (ys>y(ii+1))&(ys<y(ii));[yyy,iii]=sort([yyy,ys]);wyyy=[wyyy,0];wyyy=wyyy(iii);end; % on ajoute la source
yy=[yy,yyy];wyy=[wyy,wyyy];
[ee,oo]=deal(zeros(length(yyy),1,3));
prv1=exp(i*khi(ii).'*(yyy-y(ii+1)));prv2=exp(-i*khi(ii).'*(yyy-y(ii)));
ee(:,1,1)=a(ii,1)*prv1+a(ii,2)*prv2;
ee(:,1,2)=(i*khi(ii)/mu(ii))*(a(ii,1)*prv1-a(ii,2)*prv2);
ee(:,1,3)=(-i*g/mu(ii))*(a(ii,1)*prv1+a(ii,2)*prv2);
oo(:,1,1)=ep(ii);
oo(:,1,2:3)=mu(ii);
e=[e;ee];
o=[o;oo];
end;
[prv,a]=retreal(e(:,:,1));e=a*e;if wyy*e(:,:,1)<0;e=-e;end; % E est 'presque' réel positif

% normalisation à Lorentz=-4i (Pt=1 pour les modes propagatifs)
khih=retsqrt(ep(1)*mu(1)-g.^2,-1);khib=retsqrt(ep(end)*mu(end)-g.^2,-1);
pt=.5i*e(:,:,1).*e(:,:,3);
pt=wyy*pt-pt(1)/(2i*khib)-pt(end)/(2i*khih);
e=e/sqrt(pt);% pt=1
% calcul de vmc;
vmc=o(:,:,1).*e(:,:,1).^2+o(:,:,3).*e(:,:,2).^2-o(:,:,2).*e(:,:,3).^2;
vmc=wyy*vmc-vmc(1)/(2i*khib)-vmc(end)/(2i*khih);

e(:,:,2:end)=-i*e(:,:,2:end);    % declonage	
end;         %  <****


%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   version 1 D balayage   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [zz,kmax,z]=balayage_1D(pol,n,h,betr,beti,k0,hb);
% [zz,kmax,z]=retmode(pol,n,h,betr,beti,k0,hb);
% pol,n,h,: comme retmode 1D
% betr,beti vecteurs des neff complexes (K//=(Betr+i*Beti)*k0, où [Beti,Betr]=ndgrid(beti,betr);
% hb si on veut le tracé des modes
% Le programme calcule t avec retabeles vectorialisé.On affiche une quantité deduite de t par utilisation du theoreme de Cauchy, dont les maximas correspondent aux poles.
% On peut ensuite selectionner des points de depart z
% si zz est en sortie,les modes sont cherchés (valeur zz et kmax nombre d'iterations)

if nargin<7;hb=[];end;
n=n(:);h=h(:);
tab=[h,n(2:end-1)];
[Beti,Betr]=ndgrid(beti,betr);Beta=(Betr+i*Beti)*k0;

[rh,t,rb,prv]=retabeles(pol,n,h,Beta(:).',k0,[],0);
t=reshape(rh+t+rb,size(Beta));

S=(t(1:end-1,1:end-1)+t(1:end-1,2:end)-t(2:end,1:end-1)-t(2:end,2:end))*retdiag(diff(betr))...
-1i*(retdiag(diff(beti))*(t(1:end-1,1:end-1)+t(2:end,1:end-1)-t(1:end-1,2:end)-t(2:end,2:end)));
fig=figure;font=retfont;font{2}='tex';
X=(betr(1:end-1)+betr(2:end))/2;Y=(beti(1:end-1)+beti(2:end))/2;
X0=mean(X);Y0=mean(Y);recentre=min(std(X)/X0,std(Y)/Y0)<1.e-5;
if recentre;X0=str2num(num2str(X0,1));Y0=str2num(num2str(Y0,1));else;X0=0;Y0=0;end;
retcolor(X-X0,Y-Y0,log10(abs(S)),-2);
if recentre;set(gca,'XTickLabel',num2str(str2num(get(gca,'XTickLabel'))+X0,2),'XTickMode','manual','YTickLabel',num2str(str2num(get(gca,'YTickLabel'))+Y0,2),'YTickMode','manual');end;
hold on;
teta=retsqrt;
% cp1=sqrt(n(1)^2-linspace(0,2*(max(abs(betr))^2+max(abs(beti))^2),10000)*exp(-1i*teta));f=find(real(cp1)<=max(betr)&real(cp1)>=min(betr)&imag(cp1)<=max(beti)&imag(cp1)>=min(beti));plot(real(cp1(f)),imag(cp1(f)),'-w','linewidth',2);f=find(real(-cp1)<=max(betr)&real(-cp1)>=min(betr)&imag(-cp1)<=max(beti)&imag(-cp1)>=min(beti));plot(real(-cp1(f)),imag(-cp1(f)),'-w','linewidth',2);
% cp2=sqrt(n(end)^2-linspace(0,2*(max(abs(betr))^2+max(abs(beti))^2),10000)*exp(-1i*teta));f=find(real(cp2)<=max(betr)&real(cp2)>=min(betr)&imag(cp2)<=max(beti)&imag(cp2)>=min(beti));plot(real(cp2(f)),imag(cp2(f)),'-w','linewidth',2);f=find(real(-cp2)<=max(betr)&real(-cp2)>=min(betr)&imag(-cp2)<=max(beti)&imag(-cp2)>=min(beti));plot(real(-cp2(f)),imag(-cp2(f)),'-w','linewidth',2);
if prod(betr([1,end]))<0;plot([-X0,-X0],beti([1,end])-Y0,'--w','linewidth',2);end;
if prod(beti([1,end]))<0;plot(betr([1,end])-X0,[-Y0,-Y0],'--w','linewidth',2);end;
xlabel('\Ree(n_e_f_f )',font{:});ylabel('\Imm(n_e_f_f )',font{:});set(gca,font{3:end});drawnow;

if nargout>0;% entree et recherche de poles
[z,zz,kmax]=deal([]);num=0;
while 1;
npoles=input('entree de n poles:   Entrer n puis choisir sur le graphe (0 pour sortir)  On peut zoomer avant ');
if isempty(npoles);npoles=input('');end;
npoles,if npoles==0;break;end;
figure(fig);
[x,y,prv]=ginput(npoles);

z=[z,x(:).'+1i*y(:).'+X0+1i*Y0];
for ii=1:npoles;num=num+1;plot(x(ii),y(ii),'or');text(x(ii),y(ii),[' \leftarrow ',int2str(num)],'FontSize',18);end;
end;
if nargout>=1 & ~isempty(z);% recherche et trace des poles
[zz,kmax]=retmode(pol,n,h,z*k0,[],hb,[],k0);zz=zz/k0;zz(kmax>=100)=nan+1i*nan;
figure(fig);plot(real(zz)-X0,imag(zz)-Y0,'*c');
end;             % recherche et trace des poles
end;% entree et recherche de poles


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   version 1 D suivi de mode   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [lld,ggama,coupure]=suit_mode(pol,indice,h,hb,ld1,gama1,parm);
% [lld,ggama,coupure]=suit_mode(pol,indice,h,hb,ld1,gama1,parm);
% indice:fonction	 n=indice(ld) est le parametre de retmode
% pol,h, parametres de retmode
% hb
% ld1,gama1=K//1=neff0*k0 point de départ approximatif
% le mode de depart est cherché et eventuellement tracé (selon hb)
% 
% parm=struct('gama_reel',0,'antizap',1.5,'tracer',0, ...
% 	'ldmin',-inf,'ldmax',inf,'dld0',1.e-6,'dldmax',5.e-3,'dldmin',1.e-11, ...
% 	'gamamin',-inf,'gamamax',inf,'dgama0',1.e-6,'dgamamax',5.e-3,'dgamamin',1.e-11);
%                       par defaut
% gama_reel=0 on pole en gama=K a ld reel fixé Indice est une fonction complexe de variable reelle
% gama_reel=1 on pole en ld a gama reel fixé Indice est une fonction complexe de variable complexe (modele de type Drude necessaire)
% ldmin, ldmax bornes en ld ( ou gamamin, gamamax bornes en gama)
% dld0 (ou dgama0)pas initial (ou dgama0)
% dldmax (ou dgamamax)  plus grand pas permis
% dldmin (ou dgamamin)  plus petit pas permis ( coupure)
% Le mode est suivi dans les 2 sens depuis ld1 (ou gama1)
% au cours du calcul, le pas peut être augmenté (si kmax<5) jusqu'à dldmax (ou dgamamax)
% Quand on s'approche des coupures de mode, le pas diminue jusqu'à dldmin (ou dgamamin)
% Si la variation est trop brusque ( parametre antizap>1) , on diminue le pas
%  en sortie:
% lld longueurs d'ondes  ggama= K//=neff*k0
% coupure vecteur de 2 elements: si le mode arrive à un cut_off à gauche coupure(1)=1, si le mode arrive à un cut_off à droite coupure(2)=1
% si indice donne un cellarray, c'est en fait epsilon et on utilise retyeh

pparm=struct('gama_reel',0,'antizap',1.5,'tracer',0, ...
	'ldmin',-inf,'ldmax',inf,'dld0',1.e-6,'dldmax',5.e-3,'dldmin',1.e-11, ...
	'gamamin',-inf,'gamamax',inf,'dgama0',1.e-6,'dgamamax',5.e-3,'dgamamin',1.e-11);
gama_reel=retoptimget(parm,'gama_reel',pparm);
antizap=retoptimget(parm,'antizap',pparm);
tracer=retoptimget(parm,'tracer',pparm);

if gama_reel==1;% lld complexe, gama reel
clear fonc_mode_bis;
gamamin=retoptimget(parm,'gamamin',pparm);
gamamax=retoptimget(parm,'gamamax',pparm);
dgama0=retoptimget(parm,'dgama0',pparm);
dgamamax=retoptimget(parm,'dgamamax',pparm);
dgamamin=retoptimget(parm,'dgamamin',pparm);
[ggama,lld,coupure]=ret_suit_mode(gama1,ld1,@fonc_mode_bis,struct('xmin',gamamin,'xmax',gamamax,'dx0',dgama0,'dxmax',dgamamax,'dxmin',dgamamin,'kmax0',5,'antizap',antizap,'tracer',tracer),indice,pol,h,hb);
else;           % lld reel, gama complexe
clear fonc_mode;
ldmin=retoptimget(parm,'ldmin',pparm);
ldmax=retoptimget(parm,'ldmax',pparm);
dld0=retoptimget(parm,'dld0',pparm);
dldmax=retoptimget(parm,'dldmax',pparm);
dldmin=retoptimget(parm,'dldmin',pparm);
	
[lld,ggama,coupure]=ret_suit_mode(ld1,gama1,@fonc_mode,struct('xmin',ldmin,'xmax',ldmax,'dx0',dld0,'dxmax',dldmax,'dxmin',dldmin,'kmax0',5,'antizap',antizap,'tracer',tracer),indice,pol,h,hb);
end;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [gama,kmax,test,A]=fonc_mode(ld,gama1,indice,pol,h,hb);
persistent nv;
if ~isempty(nv);hb=[];else;nv=1;end;% on trace le champ au premier passage
k0=2*pi/ld;n=indice(ld);
if iscell(n);% yeh
for kkk=0:3;
[gama,kmax,prv,prv,test]=retcadilhac(@calg,struct('nout',4,'niter',100,'nitermin',0,'tol',1.e-10,'tolf',1.e-5),gama1*[1,1.001,.999],n,k0,h,pol+1i*kkk);
if all(test);break;end;
end;
A=1;
else;
[gama,kmax]=retmode(pol,n,h,gama1,[],hb,[],k0);test=kmax<100;
if iscell(h);A=retsqrt((k0*[n{1}(1);n{end}(end)]).^2-gama^2,-1);else;A=retsqrt((k0*[n(1);n(end)]).^2-gama^2,-1);end;
end;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [ld,kmax,test,A]=fonc_mode_bis(gama1,ld,indice,pol,h,hb);
persistent nv;
if ~isempty(nv);hb=[];else;nv=1;end;% on trace le champ au premier passage
k0=2*pi/ld;
for kkk=0:3;
[k0z,kmax,prv,prv,test]=retcadilhac(@calg_bis,struct('niter',100,'nitermin',0,'tol',1.e-10,'tolf',1.e-5),k0*[1,1.001,.999],indice,gama1,h,pol+1i*kkk);
if all(test);break;end;
end;
ld=2*pi/k0z;n=indice(ld);test=all(test);
if ~isempty(hb) & ~iscell(n) ;[gama,kmax]=retmode(pol,n,h,gama1,[],hb,[],k0z);end;
if iscell(n);A=retsqrt((k0*[n{1}(1);n{end}(end)]).^2-gama1^2,-1);else;A=retsqrt((k0*[n(1);n(end)]).^2-gama1^2,-1);end;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function gg=calg_bis(k0,indice,g,h,pol)
ld=2*pi/k0;n=indice(ld);
if iscell(n);% yeh n, est eps
gg=calg_yeh(k0,n,g,h,pol);
else;
[prv,prv,prv,gg]=calg(g,n,k0,h,pol);
end;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function gg=calg_yeh(k0,n,g,h,pol);
s=retyeh(n,h,[0,g],k0);
switch pol;
case 0;gg=1/s(1,3);
case 1i;gg=1/s(1,1);
case 2i;gg=1/s(3,1);
case 3i;gg=1/s(3,3);
	
case 2;gg=1/s(2,4);
case 2+1i;gg=1/s(2,2);
case 2+2i;gg=1/s(4,2);
case 2+3i;gg=1/s(4,4);
end;	

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   version 1 D avec modes de Bloch aux extrémités   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [gama,kmax,vm,erm,e,o,yy,wyy,vmc]=ret1mode_bloch(calvm,pol,n,h,ga,gb,h_b,ys,k0);
vm=[];erm=[];ept=[];vmc=[];
nh=n{1};hh=h{1};
nc=n{2};hc=h{2};
nb=n{3};hb=h{3};
if nargin<5;ga=[];end;if nargin<6;gb=[];end;if nargin<7;h_b=[];end;if nargin<8;ys=0;end;if nargin<9;k0=1;end;
if isempty(ys);ys=0;end;
tracer=(nargin>=7)&(~isempty(h_b));if (length(h_b)==3)&(h_b(3)==0) tracer=0;end;
if isempty(h_b);h_b=[0,0];end;
if sum(hh)==0;hh=h_b(1);if hh==0;hh=.1;end;end;if sum(hb)==0;hb=h_b(2);if hb==0;hb=.1;end;end;% pour les limites non periodiques
[gama,vmc]=deal(zeros(size(ga)));kmax=gama;[vm,erm,ept]=deal(zeros(3,length(ga)));[e,o,yy,wyy]=deal(cell(1,length(ga)));
tabh=[hh(:),nh(:)];tabc=[hc(:),nc(:)];tabb=[hb(:),nb(:)];
if isempty(ga); % trace de rh rb t
figure;
if length(gb)==3;nn=gb(3);else nn=1000;end;
if isempty(gb);x=linspace(min([real(nh*k0),real(nb*k0)]),max(real(nc*k0)),nn);else;x=linspace(gb(1),gb(2),nn);end;
init={pol,x,k0};
if size(tabb,1)==1;sb=retb(init,tabb(1,2),-1);else;[sh,sb]=retbloch(init,retcouche(init,tabb),sum(real(tabb(:,1))));end% sh est ensuite ecrasé attention à l'ordre ...
if size(tabh,1)==1;sh=retb(init,tabh(1,2),1);else;sh=retbloch(init,retcouche(init,tabh),sum(real(tabh(:,1))));end;

s=reteval(retss(sh,retcouche(init,tabc),sb));
rh=squeeze(s(1,2,1,:));rb=squeeze(s(2,1,1,:));t=squeeze(s(1,1,1,:));
plot(x,abs(rh).^2,'-k',x,abs(rb).^2,'--k',x,abs(t).^2,':k','linewidth',2);legend('rh','rb','t');xlabel(['g    deltag=',num2str(x(end)-x(1)),' g= ',num2str(.5*(x(1)+x(end)),10)]);

else % recherche des modes 
[gama,kmax]=deal(zeros(size(ga)));
for ii=1:length(ga);
[gama(ii),kmax(ii)]=retcadilhac(@calg_bloch,struct('nout',4,'niter',100,'nitermin',4,'tol',1.e-10,'tolf',1.e-5),ga(ii),nc,k0,hc,tabh,tabb,pol);
%[gama(ii),kmax(ii)]=retcadilhac(@calg_bloch,struct('nout',4,'niter',100,'nitermin',4,'tol',1.e-10,'tolf',1.e-5),ga(ii)^2,nc,k0,hc,tabh,tabb,pol);gama(ii)=sqrt(gama(ii));
end
end;

% tracé du champ des modes et calcul du volume modal
if tracer;figure;end;
if calvm|tracer;

kkk=ceil(sqrt(length(ga)));kkkk=0;
for ii=1:length(ga);%%%%%%%%%%%%%%
[rh,rb,t,gg,e{ii},yy{ii},wyy{ii},o{ii},Ath,Atb,Hh,Hb,H,yys,vmc(ii)]=calg_bloch(gama(ii),nc,k0,hc,tabh,tabb,pol,h_b,ys);
 Ath=abs(Ath);Atb=abs(Atb);
if calvm;xx=0;wxx=1;
[prv,erm(:,ii),ec,ener,enere,enerh,poynting]=retvm(e{ii},o{ii},{xx,yy{ii}},{wxx,wyy{ii}},[0,yys]);vm(:,ii)=prv(:);
% termes correctifs dus au prolongement exponentiel décroissant à l'infini
fy=find(yy{ii}<Hb);[prv,prv1,ec,ener,enere_b,enerh_b]=retvm(e{ii}(fy,:,:),o{ii}(fy,:,:),{xx,yy{ii}(fy)},{wxx,wyy{ii}(fy)},[0,yys]);
fy=find(yy{ii}>H-Hh);[prv,prv1,ec,ener,enere_h,enerh_h]=retvm(e{ii}(fy,:,:),o{ii}(fy,:,:),{xx,yy{ii}(fy)},{wxx,wyy{ii}(fy)},[0,yys]);

eneree=enere+enere_b*Atb/(1-Atb)+enere_h*Ath/(1-Ath);
enerhh=enerh+enerh_b*Atb/(1-Atb)+enerh_h*Ath/(1-Ath);

erm(:,ii)=erm(:,ii)*sqrt((enere+enerh)/(eneree+enerhh));
vm(1,ii)=vm(1,ii)*(eneree/enere);
vm(2,ii)=vm(2,ii)*(enerhh/enerh);
vm(3,ii)=vm(3,ii)*(eneree+enerhh)/(enere+enerh);

else vm=[];erm=[];yys=nan;end;

if tracer;
kkkk=kkkk+1;retsubplot(kkk,kkk,kkkk);
if length(h_b)==2;
plot(yy{ii},real(e{ii}(:,1,1)),'-k',yy{ii},imag(e{ii}(:,1,1)),'--k',yys,0,'*k','linewidth',2);leg=legend('real','imag');
else;
plot(yy{ii},abs(e{ii}(:,1,1)).^2,'-k',yys,0,'*k','linewidth',2);leg=legend('intensite');
end;
set(leg,'fontsize',max(7,7-kkk));
xlabel([' gama=',num2str(gama(ii),10),' kmax=',num2str(kmax(ii))],'fontsize',max(7,7-kkk));
hold on;axis tight;ax=axis;oo=real(sqrt((o{ii}(:,1,1).*o{ii}(:,1,2))));oo=ax(3)+.5*oo*(ax(4)-ax(3))/max(abs(oo));plot(yy{ii},real(oo),':k','linewidth',2);hold off;axis tight;
end
end;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [rh,rb,t,gg,e,yy,wyy,o,Ath,Atb,Hh,Hb,H,ys,vmc]=calg_bloch(g,n,k0,h,tabh,tabb,pol,h_b,ys);
if nargout<5;% <****   calcul de r et t
n=[tabh(end,2);n(:);tabb(1,2)];h=h(:);
[M,m,y,khi,ep,mu,ASBh,ASBb]=calM(g,n,k0,h,tabh,tabb,pol);
a=-M(:,[1:m-1,m+2:2*m])\M(:,[m,m+1]);
rh=a(1,2);
rb=a(2*m-2,1);
t=a(1,1);
gg=1/t;
gg=1/rh;
else;        %  <****   calcul du champ

nreph=ceil(h_b(1)/sum(real(tabh(:,1))));nrepb=ceil(h_b(2)/sum(real(tabb(:,1))));

ys=ys-h_b(2)+nrepb*sum(real(tabb(:,1)));% répétition des modes de Bloch	
n=[tabh(end,2);repmat(tabh(:,2),nreph,1);n(:);repmat(tabb(:,2),nrepb,1);tabb(1,2)];h=[repmat(tabh(:,1),nreph,1);h(:);repmat(tabb(:,1),nrepb,1)];

[M,m,y,khi,ep,mu,ASBh,ASBb,Dh,Db]=calM(g,n,k0,h,tabh,tabb,pol);
Hh=nreph*sum(real(tabh(:,1)));Ath=exp(2*Hh*Dh(1));% attenuation en intensite sur nreph periodes
Hb=nrepb*sum(real(tabb(:,1)));Atb=exp(2*Hb*Db(1));% attenuation en intensite sur nreph periodes
H=sum(h);

[rh,rb,t,gg,yy,wyy]=deal([]);
st_warning=warning;warning off;[a,prv]=eigs(M(:,[1:m-1,m+2:2*m]),1,'sm',struct('disp',0));
a=reshape([a(1:m-1);0;0;a(m:2*m-2)],m,2);

%y(1)=y(1)+h_b(1);y(m+1)=y(m+1)-h_b(2);
y=y-y(end);
e=zeros(0,1,3);o=zeros(0,1,3);
for ii=m:-1:1;
	
[yyy,wyyy]=retgauss(y(ii+1),y(ii),20,-1*min(100,ceil(4*abs(khi(ii))*abs(y(ii)-y(ii+1)))));
%if isempty(yyy);yyy=zeros(1,0);wyyy=zeros(1,0);end;
if (ys>y(ii+1))&(ys<y(ii));[yyy,iii]=sort([yyy,ys]);wyyy=[wyyy,0];wyyy=wyyy(iii);end; % on ajoute la source
yy=[yy,yyy];wyy=[wyy,wyyy];
[ee,oo]=deal(zeros(length(yyy),1,3));
prv1=exp(i*khi(ii).'*(yyy-y(ii+1)));prv2=exp(-i*khi(ii).'*(yyy-y(ii)));
ee(:,1,1)=a(ii,1)*prv1+a(ii,2)*prv2;
ee(:,1,2)=(i*khi(ii)/mu(ii))*(a(ii,1)*prv1-a(ii,2)*prv2);
ee(:,1,3)=(-i*g/mu(ii))*(a(ii,1)*prv1+a(ii,2)*prv2);
oo(:,1,1)=ep(ii);
oo(:,1,2:3)=mu(ii);
e=[e;ee];
o=[o;oo];
end;

[prv,a]=retreal(e(:,:,1));e=a*e;if wyy*e(:,:,1)<0;e=-e;end; % E est 'presque' réel positif
% normalisation à pt=1
pt=.5i*e(:,:,1).*e(:,:,3);
[fh,fc]=retfind(yy>H-Hh);[fb,fc]=retfind(yy<Hb,fc);
pt=wyy(1,fc)*pt(fc,1)+(wyy(1,fb)*pt(fb,1))/(1-Atb)+(wyy(1,fh)*pt(fh,1))/(1-Ath);
e=e/sqrt(pt);% pt=1
% calcul de vmc;
vmc=o(:,:,1).*e(:,:,1).^2+o(:,:,3).*e(:,:,2).^2-o(:,:,2).*e(:,:,3).^2;
vmc=wyy(1,fc)*vmc(fc,1)+(wyy(1,fb)*vmc(fb,1))/(1-Atb)+(wyy(1,fh)*vmc(fh,1))/(1-Ath);
e(:,:,2:end)=-i*e(:,:,2:end);% declonage

end;         %  <****

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [M,m,y,khi,ep,mu,ASBh,ASBb,Dh,Db]=calM(g,n,k0,h,tabh,tabb,pol);
init={pol,g,k0};
if size(tabh,1)==1;Dh=[-1;1]*retsqrt(g^2-(k0*tabh(1,2))^2,1);ASBh=eye(2);
else;[prv1,prv2,Dh,ASBh]=retbloch(init,retcouche(init,tabh),sum(real(tabh(:,1))),tabh(end,2));end;
if size(tabb,1)==1;Db=[-1;1]*retsqrt(g^2-(k0*tabb(1,2))^2,1);ASBb=eye(2);
else;[prv1,prv2,Db,ASBb]=retbloch(init,retcouche(init,tabb),sum(real(tabb(:,1))),tabb(1,2));;end;


n=n(:);h=h(:);
% graphéne
[f_pas_graphene,f_graphene]=retfind(imag(h)==0);
if ~isempty(f_graphene);
sig_G=n(f_graphene+1);
n=n([1;f_pas_graphene+1;end]);h=h(f_pas_graphene);
f_graphene=f_graphene-[0:(length(f_graphene)-1)].';% numero des dioptres à graphéne
end;


ep=retep(n,pol,k0);mu=ep(2,:).';ep=ep(1,:).';

m=length(n);
y=-[0;cumsum(h)];
khi=retbidouille(retsqrt(ep.*mu-g^2,-1));
y=[y(1);y;y(end)];
p=(1:m-1).';
prv1=exp(-i*khi(p).*(y(p+1)-y(p)));
prv2=exp(i*khi(p+1).*(y(p+1)-y(p+2)));
prv3=khi(p).*mu(p+1)./(mu(p).*khi(p+1));
[ii,jj,M]=deal([]);
ii=[ii;p];jj=[jj;p];M=[M;ones(m-1,1)];
ii=[ii;p];jj=[jj;p+m];M=[M;prv1];
ii=[ii;p];jj=[jj;p+m+1];M=[M;-ones(m-1,1)];
ii=[ii;p];jj=[jj;p+1];M=[M;-prv2];

ii=[ii;p+m-1];jj=[jj;p];M=[M;prv3];
ii=[ii;p+m-1];jj=[jj;p+m];M=[M;-prv1.*prv3];
ii=[ii;p+m-1];jj=[jj;p+m+1];M=[M;ones(m-1,1)];
ii=[ii;p+m-1];jj=[jj;p+1];M=[M;-prv2];


% graphéne
if ~isempty(f_graphene);
if pol==0;% TE
ii=[ii;f_graphene+m-1];jj=[jj;f_graphene+1];M=[M;-1i*sig_G.*prv2(f_graphene).*mu(f_graphene+1)./khi(f_graphene+1)];    
ii=[ii;f_graphene+m-1];jj=[jj;f_graphene+m+1];M=[M;-1i*sig_G.*mu(f_graphene+1)./khi(f_graphene+1)];    
else;% TM
ii=[ii;f_graphene];jj=[jj;f_graphene+1];M=[M;-1i*sig_G.*prv2(f_graphene).*khi(f_graphene+1)./mu(f_graphene+1)];    
ii=[ii;f_graphene];jj=[jj;f_graphene+m+1];M=[M;1i*sig_G.*khi(f_graphene+1)./mu(f_graphene+1)];    
end;
end;



M=sparse(ii,jj,M,2*m-2,2*m);
[M(:,[1,m+1]),M(:,[m,2*m])]=deal(M(:,[1,m+1])*ASBh,M(:,[m,2*m])*ASBb);% base des modes de Bloch 


%%%%%%%%%%%%%%%%%%%%
%   version 2 D    %
%%%%%%%%%%%%%%%%%%%%

function [gama,kmax,vm,erm,ept,e,o,xxx,zzz,wxxx,wzzz,vmc]=ret2mode(calvm,init,mm,pols,u,ga,gb,hb,fich);
if nargin<9;fich=[];end;
if isstruct(fich);parm_cadilhac=fich;else;parm_cadilhac=struct('tol',1.e-10,'tolf',1.e-5,'niter',30);end
if nargin<8;hb=[];end;
if nargin<7;gb=[];end;if isempty(gb);in=rettestobjet(init,u,-1);gb=[min(min(abs(in))),max(max(abs(in)))];end;
if nargin<6;ga=[];end;
trace=(~isempty(hb));if length(hb)>2;if(hb(3)==0);trace=0;hb=hb(1:2);end;end;

d=init{end}.d;sym=init{end}.parm.sym;cao=init{end}.parm.cao;beta0=init{end}.beta0(1);% il faut prendre parm car si pas de sym ou de cao init{end}.sym=[] ou init{end}.cao=[] 
if isempty(sym);sym=[0,0,0,0];end;sym(2)=0;
xs=sym([3,4]);

if isempty(cao);cao=[-d/2,0,0];end;cao(4)=0; %cao(3)=real(cao(3));% si on ne veut pas de cao complexe en x
if length(mm)==1;mm=[-mm,mm];end;

proteges=retio; % fichiers protégés
xs(2)=mod(xs(2),d(2));cao(2)=mod(cao(2),d(2));if cao(2)>xs(2);cao(2)=cao(2)-d(2);end;
cao(1)=mod(cao(1),d(1));if cao(1)>xs(1);cao(1)=cao(1)-d(1);end;
[wb,discb]=retautomatique([d(1),1],[],u,[d(2),cao(2),xs(2)]);wgb=wb(end);if length(wb)>1;wb=wb(1:end-1);else;wb{1}{1}=0;end;
[wh,disch]=retautomatique([d(1),1],[],u,[d(2),xs(2),cao(2)+d(2)]);wgh=wh(1);if length(wh)>1;wh=wh(2:end);else;wh{1}{1}=0;end;
wgs=wb(1);wgs{1}{1}=0;wgb{1}{1}=0;wgh{1}{1}=0;
if calvm&isempty(hb);hb=[retauto(wh),retauto(wb),0];end; % valeurs par défaut pour calcul du volume modal
[s,a,tab,lex]=retauto([],{wgs,wgh,wgb,wh,wb},-inf);


[e,o,xxx,zzz,wxxx,wzzz,vmc]=deal(cell(size(ga)));
[gama,vmc,kmax]=deal(zeros(size(ga)));vm=zeros(3,length(ga));[erm,ept]=deal(zeros(6,length(ga)));

if init{end}.granet==1;cao_granet=init{13};cao_granet(5:8)={[],[],[],[1]};cao={cao,cao_granet};end



if isempty(ga); % ********* tracé de la réponse des modes
figure;    
if length(gb)==3;nn=gb(3);else;nn=20;end
x=linspace(gb(1),gb(2),nn);r=zeros(length(pols),nn);
nf=ceil(sqrt(length(pols)));
for ii=1:nn;g=x(ii);[prv1,prv2,prv3,prv4,r(:,ii)]=calq(g,d,mm,sym,beta0,cao,lex,xs(1),pols);end
for jj=1:length(pols);
retsubplot(nf,nf,jj);plot(x,abs(r(jj,:)).^2,'-k','linewidth',2);xlim([min(x),max(x)]);
xlabel(['pols=',num2str(pols(jj)),' deltag=',num2str(x(end)-x(1)),' g= ',num2str(.5*(x(1)+x(end)),10)]);
retfont(gcf,0);drawnow;
end;

else;           %******** recherche des modes 
for iii=1:length(ga);gama(iii)=ga(iii);% iii
	
    if ischar(fich);km=30;		% ancienne version obsolete
	er=inf;k=0;z=[];zz=[];
	while(er>1.e-10)&(k<km);
	zz0=calq(gama(iii),d,mm,sym,beta0,cao,lex,xs(1),pols);
	k=k+1;z=[z;gama(iii)];zz=[zz;zz0];
	[gama(iii),z,zz,er]=retcadilhac(z,zz);
	retio(proteges,-4);
	save(fich,'gama','z','zz','er','k');
	end;
else;
    [gama(iii),k,erz,erfonc,test]=retcadilhac(@calq,parm_cadilhac,gama(iii)*(1+[-1.e-6;0;1.e-6]),d,mm,sym,beta0,cao,lex,xs(1),pols,proteges);
end;

if ~all(test);gama(iii)=0;% mode non trouvé
else;                     % calcul du champ
if (~isempty(hb));
[prv,vm(:,iii),erm(:,iii),ept(:,iii),qq,e{iii},o{iii},xxx{iii},zzz{iii},wxxx{iii},wzzz{iii},vmc(iii)]=calq(gama(iii)+1.e-8*real(gama(iii))+i*1.e-8*abs(imag(gama(iii))),d,mm,sym,beta0,cao,lex,xs(1),pols,proteges,iii*trace,calvm,hb,[disch{1},discb{1}],[' gama=',num2str(gama(iii))]);retio(proteges,-4);
end;
end;
kmax(iii)=k;
end;                                    % iii
end;      %********
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [q,vm,erm,ept,qq,e,o,x,z,wx,wz,vmc]=calq(g,d,mm,sym,beta0,cao,lex,xs,pols,proteges,trace,calvm,hb,disc,texte);
if nargin<11;trace=0;calvm=0;end;% le beta0 introduit le 8 5 2015
apod=1;7.25;
q=0;vm=0;erm=0;ept=0;
init=retinit([d(1),1],[mm,0,0],[beta0,g],sym,cao);
%                  [s,a,tab,lex]=retauto([],{wgs,wgh,wgb,wh,wb},-inf);
if (trace==0)&(~calvm);[s,a,tab]=retauto(init,[],{[1+i,xs,0],0},lex);
else;[s,a,tab]=retauto(init,[],{1+i,1},lex);end;
ks=tab{1}(2);as=a{ks};kgh=tab{2}(2);ah=a{kgh};kgb=tab{3}(2);ab=a{kgb};ssh=s{4};ssb=s{5};
sh=retb(init,ah,1.e-1,0);sh=rettronc(sh,[],[],1);
sb=retb(init,ab,-1.e-1,0);sb=rettronc(sb,[],[],-1);
source=rets(init,[xs,0],as,pols,pols,[apod,apod]);
s=retss(sh,ssh,source,ssb,sb);
qq=diag(reteval(s));
q=1/sum(qq);
if (trace>0)|calvm; % <<<<<<<<<<trace ou calcul de vm
vmc=0;
if iscell(cao);cao=cao{1};end;% granet
tabh=tab{4};tabb=tab{5};a=[a,{source}];ksource=length(a);
[x,wx]=retgauss(cao(1),d(1)+cao(1),15,10,disc,d(1));
if calvm; % calcul en plus de points
[x,wx]=retgauss(cao(1),d(1)+cao(1),15,40,disc,d(1));
tab=[[hb(1),kgh,0];tabh;[0,ks,i*ksource];tabb;[hb(2),kgb,0]];f=find(imag(tab(:,3))==0);tab(f,3)=ceil(400*tab(f,1)/sum(tab(f,1)));
[x,prv]=sort([x,xs]);wx=[wx,0];wx=wx(prv);% on rajoute la source
else
[x,wx]=retgauss(cao(1),d(1)+cao(1),15,10,disc,d(1));
tab=[[hb(1),kgh,0];tabh;[0,ks,i*ksource];tabb;[hb(2),kgb,0]];f=find(imag(tab(:,3))==0);tab(f,3)=ceil(100*tab(f,1)/sum(tab(f,1)));
%tabh(:,3)=10;tabb(:,3)=10;tab=[[hb(1),kgh,5];tabh;[0,ks,i*ksource];tabb;[hb(2),kgb,5]];
end;
y=0;wy=1;[e,z,wz,o]=retchamp(init,a,sh,sb,ones(size(pols)),{x,y},tab,[],[1:6]+i*apod,1,1,1:6);
e(:,:,:,4:6)=i*e(:,:,:,4:6);% clonage

[prv,a]=retreal(e(:,:,:,[1,3,5]));e=a*e;if wz*(e(:,:,:,1)+e(:,:,:,3)+e(:,:,:,5))*wx.'<0;e=-e;end;% Ez Ex Hy sont 'presques' réels positifs
pt=-.5i*wz*squeeze(e(:,:,:,3).*e(:,:,:,4)-e(:,:,:,1).*e(:,:,:,6))*wx.';
e=e/sqrt(pt);% pt=1
% calcul de vmc
vmc=wz*(o(:,:,:,4).*e(:,:,:,1).^2-o(:,:,:,5).*e(:,:,:,2).^2+o(:,:,:,6).*e(:,:,:,3).^2 ...
-o(:,:,:,1).*e(:,:,:,4).^2+o(:,:,:,2).*e(:,:,:,5).^2-o(:,:,:,3).*e(:,:,:,6).^2)*wx.';

% vmc=2*wz*(-o(:,:,:,1).*e(:,:,:,4).^2+o(:,:,:,2).*e(:,:,:,5).^2-o(:,:,:,3).*e(:,:,:,6).^2)*wx.';


e(:,:,:,4:6)=-i*e(:,:,:,4:6);% declonage



if (trace~=0);oo=o;[prv,nx]=min(abs(x-xs));nz=find(imag(tab(:,3))>0);nz=sum(tab(nz+1:end,3));oo(max(1,nz-1):min(nz+1,end),max(1,nx-1):min(nx+1,end),1,:)=nan;rettchamp(e,oo,x,y,z,[1:6,1i],[],[],texte);end;
if  calvm;  % calcul du volume modal
is=find(imag(tab(:,3))~=0);zs=sum(tab(is:end,1),1);
[vm,erm,ec,ener,enere,enerh,poynting]=retvm(e,o,{x,y,z},{wx,wy,wz},[xs,0,zs]);vm=vm(:);erm=erm(:);
erm=erm([1,3,2,4,6,5]);erm([3,6])=-erm([3,6]);% car les axes ne sont pas les memes
ept=erm*sqrt(enere+enerh)/sqrt(-poynting(2,1));
else
vm=zeros(3,1);erm=zeros(6,1);
end;
else;
if nargin>9;retio(proteges,-4);end;
end;    % <<<<<<<<<<trace ou calcul de vm

