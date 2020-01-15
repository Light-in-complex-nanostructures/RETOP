function varargout=retabeles(varargin);
%  function [rh,th,rb,tb]=retabeles(pol,n,h,beta,k0,teta,parm);
%  Le seul programme 'couches_ minces' qui n'utilise pas les matrices d'Abeles
%
% pol 0 E//   2 H//
% n indices de haut en bas
% h hauteurs de haut en bas (length(h)=length(n)-2
%
% si parm=1 (ou absent):
% beta: tableau des beta=nh*sin(angle d'incidence) (K//=k0*beta)
%    (khi=sqrt(k0^2*(n^2-beta^2)  avec la determination teta   ) 
% si parm=0 :
% beta: tableau des  K// 
%    (khi=sqrt((k0*n)^2-beta^2)  avec la determination teta   ) 
%
% k0: tableau des k0=2*pi/lambda ( par defaut k0=1)
%  teta determination racine carrée (facultatif par defaut la coupure actuelle (obtenue par retsqrt) )
%
% calcule les coefficients de reflexion et transmission en haut (rh th) et en bas (rb tb) 
%  d'un empilement de couches homogenes d'indice n et de hauteur h
%  rh,tb,rb,de dimension: [length(k0(:)),length(beta(:))]
%    si on  demande  rh th rb ET PAS tb, th est normalise .La transmission en energie est  T=abs(th).^2  (attention il faut aussi demander rb)
%    si on  demande  tb  ET  th  ce sont les facteurs de transmission 
%    d' ondes dont le champ (E en E//  H en H // ) a l'origine est 1. 
%    La transmission en energie est alors T=abs(th.*tb)   ('double pesee')
%
%%%  Exemple
%   incidence=linspace(-90,90,180);[rh,th,rb,tb]=retabeles(2,[1,3,1.5],[1],sin(incidence*pi/180),2*pi/1.5);
%   figure;plot(incidence,abs(rh).^2,incidence,abs(th.*tb),'--',incidence,abs(th.*tb)+abs(rh).^2,':');
%   legend('R','T','R+T');xlabel('angle d''incidence');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    %    O D     VECTORIALISE         %
%    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
%  beta: tableau des  K// =k0*nh * sin(angle d' incidence)   vecteur de longueur nb
%    (khi=sqrt((k0*n)^2-beta^2)  avec la determination teta   ) 
% k0=2*pi/ld vecteur de longueur nk
% pol: 0 E//  2 H//
% 
%  init={pol,beta,k0};  reste le mëme dans tout le problème
% s=retcouche(init,tab);   construction d'une matrice de passage dans un tronçon
% tab=[[h1,n1];[h2,n2];[h3,n3];......]; hauteurs et indices de haut en bas 
%      haut                     bas 
% ( s est un cell array de dimension 2,2 contenant les tableaux de dimension nk,nb
% correspondant à une matrice de passage entre champs E et: i *k0 * Hx = 1/n^pol d(E)/dy
% En pratique on n'aura pas à accéder à ces matrices.
% Toutes les opérations sur ces matrices sont faites vectoriellement  )
% Conditions aux limites:
%   en haut:sh=retb(init,nh,1);  passage champs  --> ondes planes 
%   en bas: sb=retb(init,nb,-1); passage  ondes planes --> champs
%    ( Les ondes planes propagatives ont un flux du vecteur de poynting egal à .5)
%    si nh ou nb=inf metal electrique (E=0 si pol=0 H=0 si pol=2) si nh ou nb=-inf metal magnetique (H=0 si pol=0 E=0 si pol=2)
% Modes de Bloch:
%    [sh_bloch,sb_bloch,d_bloch,d,ASB]=retbloch(init,s,h,n_bloch); 
%   s matrice de passage du tronçon de hauteur h ,n_bloch:indice du milieu ou on veut calculer ASB (falcultatif)
%    sh_bloch :passage champs  --> modes de bloch  
%    sb_bloch :passage  modes de bloch --> champs 
%     ( Les modes de bloch propagatifs ont un flux du vecteur de poynting egal à .5)
%    d(2,nk,nb)   valeurs propres des modes du milieu equivalent (a 2*pi*i/h prés..)
%      (  d(1,:,:) : vers le haut,d(2,:,:) : vers le bas )
%    ASB(2,2,nk,nb),composantes sur les ondes planes (vers le haut ,vers le bas )du milieu d'indice n_bloch 
%           des 2 modes de Bloch ( vers le haut ,vers le bas ) 
%  
% produits de matrices avec retss:  ss=retss(sh,s,sb_bloch)
% Pour acceder au resultat: SS=reteval(ss) est une matrice s de dimension [2,2,nk,nb]
% | Dh |                 | Ib |
% |    | = S(:,:,ik,jb) *|    |
% | Db |                 | Ih |
% Cette matrice S peut se transformer en matrice T : T=retgs(S,2)
% 
% 
% Possibilite d'inverser un tronçon de matrice s par s_inv=retrenverse(s)  ( 4 cas possibles voir retrenverse)
% Possibilite d'elever s à une puissance entière ( même <=0) avec retsp(s)
%
%  CALCUL DES CHAMPS 
% [e,y,w,o]=retchamp(init,tab,sh,sb,inc,npts);
% 
%  y vecteur de longueur ny des points où le champ est calculé 
%     par défaut le nombre de points est géré par le programme et peut être multiplié par le facteur facultatif npts
%     On peut imposer les valeurs de y en mettant npts={y}
%  w poids pour une integration (gauss) 
%  e(ny,1,3, nk , nb ) la seconde dimension 1 assure la compatibilité avec le 1 D 
%  e(iy,1, 1 ,ik,ib)  = E si  pol=0    H si pol=2
%  e(iy,1, 2 ,ik,ib)  = Hx si pol=0  -Ex si pol=2
%  e(iy,1, 3 ,ik,ib)  = Hy si pol=0  -Hy si pol=2
%  On peut tracer le champ avec:   x=0; rettchamp(e(:,:,:,ik,ib),o,  x  ,y,pol); 
% 		%
% 		% modif 12 2011 variation des indices fonction de k0 OBSOLETE
% 		% global retabeles_change_indices; où la fonction retabeles_change_indices
% 		% n_de_k0=retabeles_change_indices(k0,n_bidon);
% 		% retabeles_change_indices=@change_indices;
% 		% % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 		% % function n_nv=change_indices(k0,n);
% 		% % ld=2*pi./k0(:);
% 		% % n_nv=repmat(n(:).',length(k0),1);
% 		% % f=find(n==4.610);if ~isempty(f);n_nv(:,f)=repmat(real(retindice(ld,4.6,.1)),1,length(f));end;% Al0.1Ga0.9As
% 		% % f=find(n==4.7);if ~isempty(f);n_nv(:,f)=repmat(real(retindice(ld,4.7,'linear')),1,length(f));end;    % GaAs
% modif 8 2012 variation des indices fonction de k0
%  retabeles(num,fonc1,fonc2,...);
% num tableau de numeros d'indices 'bidons'
% fonc1(ld) --> indice à remplacer pour num(1)
% fonc2(ld) --> indice à remplacer pour num(2)
% ...		
%
% Pour effacer la variable globale  retabeles_change_indices, faire retabeles sans argument
%
% % EXEMPLES 
% n2=1;n3=3.48;h2=.183502;h3=.35-h2;h=h2+h3;
% k0=linspace(2,15,20000);beta=0;pol=0;init={pol,beta,k0};
% sh=retb(init,n3,1); sb=retb(init,n3,-1);         % calcul de sh,sb pour l'indice n3
% tab=[[4*h2/5,n2];[h3,n3];[h2/5,n2];];	s=retcouche(init,tab);   % une couche n2,n3
% %  construction progressive d'un miroir de bragg
% fig=figure;hold on;axis([min(k0),max(k0),0,1.1]);xlabel('k0');title('REFLEXION D''UN MIROIR DE BRAGG');
% for ii=[0:5,128,1024];S=reteval(retss(sh,retsp(s,ii),sb));plot(k0,squeeze(abs(S(1,2,:,:))).^2,'-g');pause(.1);end
% % etude du miroir et des modes de Bloch
% [sh_bloch,sb_bloch,D,ASB]=retbloch(init,s,h,n3);
% S=reteval(retss(sh,s,sb_bloch));plot(k0,squeeze(abs(S(1,2,:,:))).^2,'-r');      % reflexion par dessus sur miroir de bragg
% S=reteval(retss(sh_bloch,s,sb));plot(k0,squeeze(abs(S(2,1,:,:))).^2,'-b');      % reflexion par dessous sur miroir de bragg
% [kk,uu,vg,vg3_GVD,k_vgz]=retvg(k0,cos(i*D(1,:,:)*h),h,2i);                      % diagrammes de dispersion
% 
%      % trace de ASB
% figure(fig);for ii=1:length(k_vgz);plot(real(k_vgz(ii))*[1,1],[0,1],'--k');end; % bords de gap
% figure;subplot(3,1,1);plot(k0,squeeze(abs(ASB(2,1,:)./ASB(1,1,:))),'-k',k0,squeeze(abs(ASB(1,2,:)./ASB(2,2,:))),'--r');legend('abs(ASB1)','abs(ASB2)');
% subplot(3,1,2);plot(k0,squeeze(real(ASB(2,1,:)./ASB(1,1,:))),'-k',k0,squeeze(real(ASB(1,2,:)./ASB(2,2,:))),'--r');legend('real(ASB1)','real(ASB2)');
% subplot(3,1,3);plot(k0,squeeze(imag(ASB(2,1,:)./ASB(1,1,:))),'-k',k0,squeeze(imag(ASB(1,2,:)./ASB(2,2,:))),'--r');legend('imag(ASB1)','imag(ASB2)');xlabel('k0');drawnow;
% % calcul du champ
% [e,y,w,o]=retchamp(init,[tab;tab;tab],sh_bloch,sb_bloch,[1,0]);% mode de bloch incident du bas
%   % rettchamp(e(:,:,:,1),o,0,y,pol); 
% figure('color','w');for ii=1:50:length(k0);E=abs(e(:,:,1,ii));
% subplot(2,1,1);plot(y,o/10,'-r',y,(E/max(E)).^2,'-k',[y(1),y(end)],[0,0],'linewidth',2);axis([0,max(y),0,1]);axis off;title('abs(E) ^2');
% subplot(2,2,3);plot(real(cos(i*D(1,:)*h)),k0,'-b',[-1,-1],[min(k0),max(k0)],'--g',[1,1],[min(k0),max(k0)],'--g',real(cos(i*D(1,ii)*h)),k0(ii),'or','markerfacecolor','r','markersize',6,'linewidth',2);xlabel('cos(K// h)');ylabel('k0');drawnow;
% subplot(2,2,4);plot(squeeze(abs(ASB(2,1,:)./ASB(1,1,:))),k0,'-b',abs(ASB(2,1,ii)./ASB(1,1,ii)),k0(ii),'or','markerfacecolor','r','markersize',6,'linewidth',2);xlabel('ASB');ylabel('k0');axis([0,1,-inf,inf]);drawnow;
% end;

% See also:RETSS,RETGS,RETBLOCH,RETCOUCHE,RETINVERSE,RETEVAL,RETCHAMP,RETTCHAMP,RETBLOCH,RETB

if nargin==1 & iscell(varargin{1});[varargout{1:nargout}]=retchange_indices(varargin{1}{1},varargin{1}{2}{:});return;end;
if nargin==0 | isa(varargin{2},'function_handle');[varargout{1:nargout}]=retchange_indices(varargin{:});

else;[varargout{1:nargout}]=retabeles1(varargin{:});
end;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function retchange_indices(n_bidon,varargin);
% retchange_indices(num,fonc1,fonc2,...);
% num tableau de numeros d'indices 'bidons'
% fonc1(ld) --> indice à remplacer pour num(1)
% fonc2(ld) --> indice à remplacer pour num(2)
% ...
% See also RETABELES RETINDICE RETGRAPHENE
global retabeles_change_indices;
if nargin==0;retabeles_change_indices=[];else;retabeles_change_indices=@(k0,n) retchange_indices1(k0,n,n_bidon,varargin{:});end;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function n_nv=retchange_indices1(k0,n,n_bidon,varargin);
ld=2*pi./k0(:);
n_nv=repmat(n(:).',length(k0),1);

for ii=1:length(n_bidon);
f=find(n==n_bidon(ii));if ~isempty(f);n_nv(:,f)=repmat(varargin{ii}(ld),1,length(f));end;
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [rh,th,rb,tb,T,Th,Tb,E,Y,W,O]=retabeles1(pol,n,h,beta,k0,teta,parm,nnpts);%prv={pol,n,h,beta,k0,teta,parm},pol,n,h
stabilise=1; % modif 2011
if nargin<8;nnpts=1;end;if nargin<7;parm=[];end;if nargin<6;teta=[];end;if nargin<5;k0=[];end;% par defaut
if isempty(beta);beta=0;end;if isempty(k0);k0=1;end;%if isempty(teta);teta=retsqrt;end;% modif 2 2013
if isempty(parm);parm=1;end;
rh=[];th=[];rb=[];tb=[];T=[];Th=[];Tb=[];E=[];Y=[];O=[];

n=n(:).';h=h(:).';champ=(length(n)==length(h));if champ;if isempty(n);n=[1,1];else;n=[n(1),n,n(end)];end;end;% si champ on ne calcule pas r et t 
[khi,npol]=calkhi(beta(:),k0(:),n(end),teta,pol,parm);  % en bas
v=[ones(1,length(khi));-1i*khi./npol];

log_vnorm=0;% modif 2011
if nargout>2;normb=khi./npol;w=[ones(1,length(khi));i*normb];end;
nk=length(k0(:));nb=length(beta(:));
if ((~champ)&(nargout>5)) | ((champ)&(nargout==6)); % matrice Tb
Tb=cell(2,2);	
%prv=reshape(sqrt(abs(normb)),nk,nb);Tb{1,1}=1./prv;Tb{2,1}=i*reshape(normb,nk,nb)./prv;
%prv=reshape(sqrt(khi),nk,nb)/(n(end)^(pol/2));% modif 6 2011 normalisation Lorentz
prv=reshape(sqrt(khi./npol),nk,nb);% modif 6 2011 normalisation Lorentz

Tb{1,1}=1./prv;Tb{2,1}=1i*prv;
Tb{1,2}=Tb{1,1};Tb{2,2}=-Tb{2,1};
end;

if nargout>7; % calcul du champ
if iscell(nnpts);cal_Y=0;Y=retcolonne(nnpts{1},1);Y_end=0;W=[];O=[];E=zeros(length(khi),0,2,2);
else;
cal_Y=1;E=zeros(length(khi),1,2,2);E(:,1,1,1)=1;E(:,1,2,1)=0;E(:,1,1,2)=0;E(:,1,2,2)=1;Y=0;W=0;O=n(end);
end;
end;      
%---------------------------------------------------------
global retabeles_change_indices;
for ii=length(h):-1:1;% ii *************************************
if imag(h(ii))~=0;% couche de graphene sig_G=n(ii+1)
%     if ~isempty(retabeles_change_indices);nn=retabeles_change_indices(k0,n(ii+1)).';else;nn=n(ii+1);end;% variation de sigma avec ld
% 	if pol==0;v(2,:)=v(2,:)-nn.*v(1,:).*k0;% TE a comprendre ...............
% 	if nargout>2;w(2,:)=w(2,:)-nn.*w(1,:).*k0;end;
% 	else;v(1,:)=v(1,:)+nn.*v(2,:)./k0;% TM  a comprendre ...............
% 	if nargout>2;w(1,:)=w(1,:)+nn.*w(2,:)./k0;end;
%     end;

    if ~isempty(retabeles_change_indices);nn=retabeles_change_indices(k0,n(ii+1)).';else;nn=n(ii+1);end;% variation de sigma avec ld
	if pol==0;nn_k0=retcolonne(repmat(nn(:).*k0(:),1,length(beta)),1);v(2,:)=v(2,:)-v(1,:).*nn_k0;% TE a comprendre ...............
	if nargout>2;w(2,:)=w(2,:)-w(1,:).*nn_k0;end;
	else;nnsk0=retcolonne(repmat(nn(:)./k0(:),1,length(beta)),1);v(1,:)=v(1,:)+v(2,:).*nnsk0;% TM  a comprendre ...............
	if nargout>2;w(1,:)=w(1,:)+w(2,:).*nnsk0;end;
    end;
else;
%[khi,npol]=calkhi(beta(:),k0(:),n(ii+1),teta,pol,parm); % une couche au centre 
[khi,npol]=calkhi(beta(:),k0(:),n(ii+1),0,pol,parm); % une couche au centre imag(khi)>=0 modif 26 4 2015: 0 et pas teta


vv=v;
if stabilise;% modif 2011
exp_m1=retexpm1(2i*khi*h(ii));
c=.5*exp_m1+1;s=exp_m1./(2i.*khi);f=find(abs(khi*h(ii))<10*eps);s(f)=h(ii);
else;
c=cos(khi*h(ii));s=h(ii)*retsinc(khi*h(ii)/pi);
end;	
v(1,:)=vv(1,:).*c+npol.*vv(2,:).*s;
v(2,:)=-(vv(1,:).*s.*(khi.^2))./npol+vv(2,:).*c;
if stabilise;v=v*retdiag(exp(-1i*real(khi*h(ii))));dlog_vnorm=imag(khi*h(ii));log_vnorm=log_vnorm+dlog_vnorm;end; % modif 2011
if nargout>2;
ww=w;    
w(1,:)=ww(1,:).*c+npol.*ww(2,:).*s;
w(2,:)=-(ww(1,:).*s.*(khi.^2))./npol+ww(2,:).*c;
if stabilise;w=w*retdiag(exp(-1i*real(khi*h(ii))));end; % modif 2011
end;



if nargout>7;  % calcul du champ dans la couche
%vv=v;ww=w; 
%if stabilise;vv=vv*retdiag(exp(llog_vnorm));ww=ww*retdiag(exp(llog_vnorm));vv(~isfinite(vv))=0;ww(~isfinite(ww))=0;end;% modif 2011
% test 2011;

if cal_Y; % le programme determine les Y
npts=ceil(nnpts*max(1,max(ceil((abs(real(khi))+abs(imag(khi)))*h(ii)/(.1*pi)))));
[hh,w_gauss]=retgauss(0,h(ii),npts,-1);
Y=[Y,Y(end)+hh];
else;     % les Y sont donnés
if ii>1;f=find(Y>=Y_end & Y<Y_end+h(ii));else;f=find(Y>=Y_end);end;% pour avoir le dernier point
if isempty(f);hh=[];else;hh=Y(f)-Y_end;w_gauss=zeros(1,length(hh));end;
Y_end=Y_end+h(ii);
end;

npts=length(hh);
if npts>0; % <<<---
W=[W,w_gauss];O=[O;repmat(n(ii+1),npts,1)];
if stabilise;% modif 2015 .......
    khih_hh=khi.'*hh;
expp=exp(1i*khih_hh+repmat(-dlog_vnorm.',1,length(hh)));
expm=exp(-1i*khih_hh+repmat(-dlog_vnorm.',1,length(hh)));
c=.5*(expp+expm);s=retdiag(1./(2i*khi))*(expp-expm);
% prv=retdiag(exp(llog_vnorm-log_vnorm));
% c=prv*c;s=prv*s;

prv=exp(-dlog_vnorm);
f=find(abs(khih_hh/pi)<.5);
if ~isempty(f);pprv=prv(:)*hh;s(f)=retsinc(khih_hh(f)/pi).*pprv(f);clear pprv;end;

E=reshape(retdiag(prv)*reshape(E,length(khi),[]),size(E));% mise à jour du champ
else; % .......
c=cos(khi.'*hh);s=retsinc(khi.'*hh/pi)*retdiag(hh);
end;

V1=retdiag(vv(1,:))*c+retdiag(npol.*vv(2,:))*s;V2=-(retdiag(vv(1,:).*(khi.^2)./npol))*s+retdiag(vv(2,:))*c;
W1=retdiag(ww(1,:))*c+retdiag(npol.*ww(2,:))*s;W2=-retdiag(ww(1,:).*(khi.^2)./npol)*s+retdiag(ww(2,:))*c;
E_prv=zeros(length(khi),npts,2,2);
E_prv(:,:,1,1)=.5*(V1+W1);
E_prv(:,:,2,1)=.5*(V2+W2);
E_prv(:,:,1,2)=retdiag(.5i./normb)*(V1-W1);
E_prv(:,:,2,2)=retdiag(.5i./normb)*(V2-W2);
if stabilise;E_prv(~isfinite(E_prv))=0;end;% modif 2011

E=cat(2,E,E_prv);
elseif stabilise;% <<<--- Ne pas oublier de renormaliser meme si pas de points dans la couche
E=reshape(retdiag(exp(-dlog_vnorm))*reshape(E,length(khi),[]),size(E));% mise à jour du champ
end;       % <<<---
end;   % fin calcul du champ dans la couche





end;% imag(h)==0
end; % boucle sur ii   *************************************
%---------------------------------------------------------

[khi,npol]=calkhi(beta(:),k0(:),n(1),teta,pol,parm);   % en haut
if (~champ)&(nargout>2);
rb=zeros(length(k0(:)),length(beta(:)));tb=zeros(length(k0(:)),length(beta(:)));    
rb(:)=-((i*khi.*w(1,:)./npol-w(2,:))./(i*khi.*v(1,:)./npol-v(2,:))).';
tb(:)=exp(log_vnorm.').*(rb(:).*v(1,:).'+w(1,:).');
end;


if ((~champ)&(nargout>4)) | (champ& ((nargout==5)|(nargout>7))); % matrice T
T=cell(2,2);
%if stabilise;prv=exp(log_vnorm);else;prv=1;end;% modif 2015
if stabilise & (nargout<=7);prv=exp(log_vnorm);else;
    prv=1;end;% modif 2015
T{1,1}=reshape(prv.*  (.5*(v(1,:)+w(1,:)))             ,nk,nb);
T{2,1}=reshape(prv.*  (.5*(v(2,:)+w(2,:)))             ,nk,nb);
T{1,2}=reshape(prv.*  (.5i*((v(1,:)-w(1,:))./normb))   ,nk,nb);
T{2,2}=reshape(prv.*  (.5i*( (v(2,:)-w(2,:)) ./normb)) ,nk,nb);
end; % matrice T

if ((~champ)&(nargout>5)) | ((champ)&(nargout==6)); % matrice Th
normh=khi./npol;
Th=cell(2,2);
%prv=reshape(sqrt(abs(normh)),nk,nb);Th{1,1}=.5*prv;Th{1,2}=-.5*i*prv./reshape(normh,nk,nb);
prv=reshape(sqrt(khi./npol),nk,nb);% modif 6 2011 normalisation Lorentz
Th{1,1}=.5*prv;Th{1,2}=-.5i./prv;
Th{2,1}=Th{1,1};Th{2,2}=-Th{1,2};
end;

if nargout>7; % calcul du champ
if cal_Y;E=cat(2,E,E(:,end,:,:));Y=[Y,Y(end)];W=[W,0];O=[O;n(1)];end;
E=reshape(E,length(k0(:)),length(beta(:)),length(Y),2,2);
end;

if ~champ;
rh=zeros(length(k0(:)),length(beta(:)));th=zeros(length(k0(:)),length(beta(:)));
v(2,:)=i*(v(2,:)./khi).*npol;
th(:)=2.*exp(-log_vnorm.')./(v(1,:)+v(2,:)).';
rh(:)=(v(1,:)-v(2,:)).'./(v(1,:)+v(2,:)).';

% if nargout==3;  % normalisation de th si on ne calcule pas tb
% normh=khi./npol;
% th(:)=th(:).*sqrt(normb./normh).';% normalisation
% modif 22_4_2015 Frederic
if nargout>2;% calcul de tb
normh=khi./npol;
if nargout==3;  % normalisation de th si on ne calcule pas tb
th(:)=th(:).*sqrt(normb./normh).';% normalisation
else;
tb(:)=th(:).*(normb./normh).';% on recalcule tb plus stable
end
end;
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [khi,npol]=calkhi(beta,k,n,teta,pol,parm);global retabeles_change_indices;
if ~isempty(retabeles_change_indices);n=retabeles_change_indices(k,n);end

if parm==1;
% [khi,prv]=ndgrid(n(:).^2-beta(:).^2,k(:).^2); % sqrt(k*(n^2-beta^2))
% khi=retsqrt1(khi.*prv,teta);
% [khi,prv]=ndgrid(n(:).^2-beta(:).^2,k(:).^2); % sqrt(k*(n^2-beta^2))
% khi=retsqrt1(khi.*prv,teta);
khi=retsqrt1(repmat((k(:).*n(:)).^2,1,length(beta))-(k(:).^2)*(beta(:).^2).',teta);
else;
%[khi,prv]=meshgrid((n(:).*k(:)).^2,beta(:).^2);   % sqrt((k*n)^2-beta^2)
[khi,prv]=ndgrid((n(:).*k(:)).^2,beta(:).^2);   % sqrt((k*n)^2-beta^2)
khi=retsqrt1(khi-prv,teta);
end;
khi(abs(khi(:))<eps)=eps; % bidouille
khi=khi(:).';

%npol=repmat(n.',size(beta,1),1);

%npol=npol(:).'.^pol;
if numel(n)==1;npol=n^pol;% 2012
else;    
%npol=repmat(n.',size(beta,1),1);
npol=repmat(n,1,size(beta,1));
npol=npol(:).'.^pol;
end;
%npol=n(:).'.^pol;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function z=retsqrt1(z,teta);if isempty(teta);z=retsqrt(z,-1);return;end;% modif 2 2013
if teta==pi;z=sqrt(z);else;z=retsqrt(exp(i*teta)*z)*exp(-i*teta/2);end;





