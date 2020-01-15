function varargout=ret_modes_Benjamin(varargin);

% Amplitude émise sur un mode d'un systeme stratifié n ,z 
%   n indices de haut en bas, z cotes de haut en bas 
% par plusieurs sources  S placées aux points M
%
% en 3D
%%%%%%%%%%%%%%
% S(6,nb_sources) en M(3,nb_sources),Attention à l'orientation de ces vecteurs
% angles=ret_modes_Benjamin([n,neff,pol],z,LL,teta,M,S,k0,parm);
% si LL=[] calcul direct pour les teta
% si LL est non vide calcul sur LLL exp(iLL Teta)
%     si LL=L1 ( une seule valeur)  LLL=L1, 
%     si LL=[L1,L2]   LLL=L1:L2, 
%     si LL=[L1,L2,L3], LLL=L1:L2:L3,
%     sinon LLL=retelimine(LL)
%      calcul de coefficients de Fourier LLL, puis synthése de Fourier en teta.
% angles structure à champs de champs f, amp, F f coefficients de Fourier, amp amplitude(teta) F intensite(teta), LLL
% si [teta,wteta]=retgauss(0,2*pi,...), l'energie sur le mode est sum(angles.F.*wteta(:)) ou 2*pi*sum(abs(f).^2)
% 
% Ajout 7 2013 pour Anthony :possibilité de calculer le champ des modes en mettant dans parm les champs:
% parm=struct('x',x,'y',y,'z',z,'xymesh',0);
% xymesh=0 (par defaut) si x et y ne sont pas 'meshés' alors size(e)= [length(z),length(x),length(y),6]
% xymesh=1 si x et y sont 'meshés' alors size(e)= [length(z),length(x)=length(y) ,1,6]
% e est alors un champ de la structure angles. Ne fonctionne que si L est non vide
%
% en 2D
%%%%%%%%%%%%%%
% S(3,nb_sources) en M(2,nb_sources),Attention à l'orientation de ces vecteurs
% angles=ret_modes_Benjamin([n,neff,pol],y,M,S,k0,parm);
% angles structure à champs de champs:   amp_p, Fp         ,amp_m,Fm
%                                  amplitude vers x>0  amplitude vers x<0
%
% parm=struct('clone',1);% par defaut
%
% See also RETCARMINATI RETBAHAR RETOP

if (nargin==6 & ~isstruct(varargin{end})) |(nargin==7 & isstruct(varargin{end}))  ;[varargout{1:nargout}]=ret_modes_Benjamin_conique(varargin{:});return;end;
if nargin>=7;[varargout{1:nargout}]=ret_modes_Benjamin_3D(varargin{:});
else;[varargout{1:nargout}]=ret_modes_Benjamin_2D(varargin{:});end;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function angles=ret_modes_Benjamin_conique(n,y,M,SS,gama_conique,k0,parm);% gama_conique du conique, gama du mode
if nargin<7;parm=[];end;
defaultopt=struct('clone',1);% par defaut la source est considéree comme clonée 
clone=retoptimget(parm,'clone',defaultopt,'fast');
if clone==0;SS(1:3,:)=1i*SS(1:3,:);end; %clonage de la source qui est non clonée
[neff,pol,n]=deal(n(end-1),n(end),n(1:end-2));
% remise eventuelle de toutes les sources à l'origine en z
for ii=1:size(M,2);SS(:,ii)=SS(:,ii)*exp(-1i*gama_conique*M(3,ii));M(3,ii)=0;end;

% elimination des points identiques
[M,k,kk]=retelimine(M.',i+100*eps);
M=M.';
S=zeros(6,length(k));
for ii=1:length(k);S(:,ii)=sum(SS(:,kk==ii),2);end;
[yy,prv,kyy]=retelimine(M(2,:));% tri important pour calmode
[gama,e]=calmode(n,y,neff,pol,k0,yy);%e=reshape(e,[],3);e=e(kyy,:);%e(:,:,1)=-e(:,:,1);% orientation
e=reshape(e(kyy,:,:),[],3);%e=e(kyy,:);
ee=expand_conique(e,M(1,:),M(2,:),gama,-gama_conique,pol,-1);amp_p=sum(sum(ee.*S.'))/(-4i);
ee=expand_conique(e,M(1,:),M(2,:),gama,-gama_conique,pol,1);amp_m=sum(sum(ee.*S.'))/(-4i);
angles=struct('amp_p',amp_p,'Fp',abs(amp_p)^2,'amp_m',amp_m,'Fm',abs(amp_m)^2,'neff',gama/k0);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function e_cartesien=expand_conique(e,X,Y,gama,gama_conique,pol,sens);
X=X(:);
gama_x=retsqrt(gama^2-gama_conique^2,-1);
cost=sens*gama_x/gama;sint=gama_conique/gama;
prv=exp(1i*gama_x*X*sens);
e_cartesien=zeros(length(X),6);
if pol==0;
e_cartesien(:,1)=(-e(:,1)*sint).*prv;
e_cartesien(:,3)=(e(:,1)*cost).*prv;
e_cartesien(:,4)=(e(:,2)*cost).*prv;
e_cartesien(:,6)=(e(:,2)*sint).*prv;
e_cartesien(:,5)=e(:,3).*prv;
else;
e_cartesien(:,1)=(e(:,2)*cost).*prv;
e_cartesien(:,3)=(e(:,2)*sint).*prv;
e_cartesien(:,2)=e(:,3).*prv;
e_cartesien(:,4)=(-e(:,1)*sint).*prv;
e_cartesien(:,6)=(e(:,1)*cost).*prv;
end;
e_cartesien=reshape(e_cartesien,size(e,1),6);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function angles=ret_modes_Benjamin_2D(n,y,M,SS,k0,parm);
if nargin<6;parm=[];end;
defaultopt=struct('clone',1);% par defaut la source est considéree comme clonée 
clone=retoptimget(parm,'clone',defaultopt,'fast');
if clone==0;SS(1,:)=1i*SS(1,:);end; %clonage de la source qui est non clonée

%rdp=randperm(size(M,2));M=M(:,rdp);SS=SS(:,rdp);
[neff,pol,n]=deal(n(end-1),n(end),n(1:end-2));
% elimination des points identiques
[M,k,kk]=retelimine(M.',1i+100*eps);
M=M.';
S=zeros(3,length(k));
for ii=1:length(k);S(:,ii)=sum(SS(:,kk==ii),2);end;
[yy,prv,kyy]=retelimine(M(2,:));% tri important pour calmode
[gama,e]=calmode(n,y,neff,pol,k0,yy);%e=reshape(e,[],3);e=e(kyy,:);%e(:,:,1)=-e(:,:,1);% orientation
% [gama,e]=calmode(n,y,neff,pol,k0,M(2,:));
% e(:,:,1)=-e(:,:,1);% orientation
e=reshape(e(kyy,:,:),[],3);%e=e(kyy,:);
ee=expand_2D(e,M(1,:),gama,-1);amp_p=sum(sum(ee.*S.'))/(-4i);
ee=expand_2D(e,M(1,:),gama,1);amp_m=sum(sum(ee.*S.'))/(-4i);
angles=struct('amp_p',amp_p,'Fp',abs(amp_p)^2,'amp_m',amp_m,'Fm',abs(amp_m)^2,'neff',gama/k0);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function e_cartesien=expand_2D(e,x,gama,sens);
e_cartesien=retdiag(exp(1i*gama*x*sens))*e;
if sens==-1;e_cartesien(:,3)=-e_cartesien(:,3);end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function angles=ret_modes_Benjamin_3D(n,z,L,teta,M,SS,k0,parm);
if nargin<8;parm=[];end;
defaultopt=struct('clone',1,'x',[],'y',[],'z',[],'xymesh',0);% par defaut la source est considéree comme clonée 
clone=retoptimget(parm,'clone',defaultopt,'fast');

if clone==0;SS(1:3,:)=1i*SS(1:3,:);end; %clonage de la source qui est non clonée
[neff,pol,n]=deal(n(end-1),n(end),n(1:end-2));

% elimination des points identiques
[M,k,kk]=retelimine(M.',i+100*eps);
M=M.';
S=zeros(6,length(k));
for ii=1:length(k);S(:,ii)=sum(SS(:,kk==ii),2);end;

[zz,prv,kzz]=retelimine(M(3,:));% tri important pour calmode
[gama,e]=calmode(n,z,neff,pol,k0,zz);e=e(kzz,:,:);e(:,:,1)=-e(:,:,1);% orientation

%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~isempty(L);% bessel
% LLL=L(1):L(2);
switch length(L);% modif 4 2011
case 1;LLL=L;
case 2;LLL=L(1):L(2);
case 3;LLL=L(1):L(2):L(3);
otherwise;LLL=retelimine(L); 
end;

ee=enroule(e,M(1,:),M(2,:),-LLL,gama,pol);
f=zeros(size(LLL));
for ii=1:length(LLL);f(ii)=-.25i*sqrt(gama/(2*pi))*exp(.25i*(2*LLL(ii)-1)*pi)*sum(sum((ee{ii}.*S.')));end;
amp=f*exp(1i*LLL.'*retcolonne(teta,1));

else;% cartesien
LLL=[];
amp=zeros(size(teta));
for ii=1:length(teta);
ee=expand_3D(e,M(1,:),M(2,:),M(3,:),pi+teta(ii),gama,pol);
amp(ii)=sum(sum(ee.*S.'));
end;
amp=-sqrt(gama/(32*pi))*exp(.25i*pi)*amp;
f=[];
end;% Bessel ou cartesien

angles=struct('teta',teta,'f',f,'amp',amp,'F',abs(amp).^2,'L',LLL,'neff',gama/k0);

% calcul du champ
if isempty(LLL);return;end;
Z=retoptimget(parm,'z',defaultopt,'fast');
if isempty(Z);return;end;
X=retoptimget(parm,'x',defaultopt,'fast');
Y=retoptimget(parm,'y',defaultopt,'fast');
xymesh=retoptimget(parm,'xymesh',defaultopt,'fast');
X=X(:);Y=Y(:);Z=Z(:);
sz=length(Z);
if xymesh==0;sx=length(X);sy=length(Y);[X,Y]=ndgrid(X,Y);else;sx=length(X);sy=1;end;
X=X(:);Y=Y(:);Z=Z(:);
[prv,Y]=ndgrid(Z,Y);[Z,X]=ndgrid(Z,X);
[zz,prv,kzz]=retelimine(Z);
[gama,e]=calmode(n,z,neff,pol,k0,zz);e=e(kzz,:,:);e(:,:,1)=-e(:,:,1);% orientation
ee=enroule(e,X,Y,LLL,gama,pol,1);
e=0;for ii=1:length(f);
Fac=(-1)^LLL(ii)*gama/2*sqrt(2*pi/gama)*exp(-.25i*(2*LLL(ii)-1)*pi);
	e=e+f(ii)*Fac*ee{ii};end;
angles.e=reshape(e,sz,sx,sy,6);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [gama,e1,e2,e3]=calmode(n,y,neff,pol,k0,y1,y2,y3);
% calcul du mode 0D et du champ normalisé pour divers cotes
y=y(:);n=n(:);hrab_min=.02*pi/real(k0);% attention hrab ne doit pas etre nul pour la normalisation du mode
hrabh=hrab_min;hrabb=hrab_min;
hrabh=max(hrabh,(y1(end)-y(1))*1.1);hrabb=max(hrabb,(y(end)-y1(1))*1.1);
if nargin>6;hrabh=max(hrabh,(y2(end)-y(1))*1.1);hrabb=max(hrabb,(y(end)-y2(1))*1.1);end;
if nargin>7;hrabh=max(hrabh,(y3(end)-y(1))*1.1);hrabb=max(hrabb,(y(end)-y3(1))*1.1);end;
[gama,kmax,vm,er,emode,o,ymode]=retmode(pol,n,-diff(y),k0*neff,[],[hrabh,hrabb,0],[],k0);
tab=[[hrabh;-diff(y);hrabb],n];
init={pol,gama,k0};
sh=retb(init,n(1),1);sb=retb(init,n(end),-1);
e1=retchamp(init,tab,sh,sb,[1,0],{ymode{1}(2:end-1)});
fac=retcolonne(e1(:,:,1:2))\retcolonne(emode{1}(2:end-1,:,1:2));% pour normalisation
%figure;plot(ymode{1}(2:end-1),real(squeeze(fac*e1(:,1,:))),ymode{1},real(squeeze(emode{1}(:,1,:))),ymode{1}(2:end-1),imag(squeeze(fac*e1(:,1,:))),'--',ymode{1},imag(squeeze(emode{1}(:,1,:))),'--')
e1=retchamp(init,tab,sh,sb,[fac,0],{y1-y(end)+hrabb});e1(:,:,2:3)=1i*e1(:,:,2:3);% clonage
if nargin>6;
e2=retchamp(init,tab,sh,sb,[fac,0],{y2-y(end)+hrabb});e2(:,:,2:3)=1i*e2(:,:,2:3);% clonage
end;
if nargin>7;
e3=retchamp(init,tab,sh,sb,[fac,0],{y3-y(end)+hrabb});e3(:,:,2:3)=1i*e3(:,:,2:3);% clonage
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function e_bessel=enroule(e,X,Y,L,gama,pol,sens);
if nargin<7;sens=0;end;
e=reshape(e,[],3);X=X(:);Y=Y(:);
e_bessel=cell(1,length(L));
[e_bessel{:}]=deal(zeros(size(e,1),6));
R=sqrt(X.^2+Y.^2);
Teta=atan2(Y,X);cosTeta=cos(Teta);sinTeta=sin(Teta);
[LL,k,kL]=retelimine([L-1,L,L+1]);
[RR,k,kR]=retelimine(R);
if sens==0;
J=retbessel('j',LL,gama*RR);
else;J=retbessel('h',LL,1,gama*RR);
end;
kL=reshape(kL,[],3);
for iL=1:length(L);
if pol==0;
e_bessel{iL}(:,1)=-.5*e(:,1).*(J(kR,kL(iL,3))+J(kR,kL(iL,1)));
e_bessel{iL}(:,2)=.5i*e(:,1).*(J(kR,kL(iL,3))-J(kR,kL(iL,1)));
e_bessel{iL}(:,4)=.5i*e(:,2.).*(J(kR,kL(iL,3))-J(kR,kL(iL,1)));
e_bessel{iL}(:,5)=.5*e(:,2).*(J(kR,kL(iL,3))+J(kR,kL(iL,1)));
e_bessel{iL}(:,6)=e(:,3).*J(kR,kL(iL,2));
else;
e_bessel{iL}(:,1)=.5i*e(:,2).*(J(kR,kL(iL,3))-J(kR,kL(iL,1)));
e_bessel{iL}(:,2)=.5*e(:,2).*(J(kR,kL(iL,3))+J(kR,kL(iL,1)));
e_bessel{iL}(:,3)=e(:,3).*J(kR,kL(iL,2));
e_bessel{iL}(:,4)=-.5*e(:,1).*(J(kR,kL(iL,3))+J(kR,kL(iL,1)));
e_bessel{iL}(:,5)=.5i*e(:,1).*(J(kR,kL(iL,3))-J(kR,kL(iL,1)));
end;
for ii=1:6;e_bessel{iL}(:,ii)=e_bessel{iL}(:,ii).*exp(1i*L(iL)*Teta);end;
[e_bessel{iL}(:,1),e_bessel{iL}(:,2)]=deal(e_bessel{iL}(:,1).*cosTeta-e_bessel{iL}(:,2).*sinTeta,e_bessel{iL}(:,1).*sinTeta+e_bessel{iL}(:,2).*cosTeta);% Ex Ey
[e_bessel{iL}(:,4),e_bessel{iL}(:,5)]=deal(e_bessel{iL}(:,4).*cosTeta-e_bessel{iL}(:,5).*sinTeta,e_bessel{iL}(:,4).*sinTeta+e_bessel{iL}(:,5).*cosTeta);% Hx Hy
e_bessel{iL}=reshape(e_bessel{iL},size(e,1),6);
end; % iL

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function e_cartesien=expand_3D(e,X,Y,z,teta,gama,pol);
e=reshape(e,[],3);X=X(:);Y=Y(:);
e_cartesien=zeros(length(z),6);
cost=cos(teta);sint=sin(teta);
prv=exp(1i*gama*(cost*X+sint*Y));
if pol==0;
e_cartesien(:,1)=(-e(:,1)*sint).*prv;
e_cartesien(:,2)=(e(:,1)*cost).*prv;
e_cartesien(:,4)=(e(:,2)*cost).*prv;
e_cartesien(:,5)=(e(:,2)*sint).*prv;
e_cartesien(:,6)=e(:,3).*prv;
else;
e_cartesien(:,1)=(e(:,2)*cost).*prv;
e_cartesien(:,2)=(e(:,2)*sint).*prv;
e_cartesien(:,3)=e(:,3).*prv;
e_cartesien(:,4)=(-e(:,1)*sint).*prv;
e_cartesien(:,5)=(e(:,1)*cost).*prv;
end;
e_cartesien=reshape(e_cartesien,size(e,1),6);







