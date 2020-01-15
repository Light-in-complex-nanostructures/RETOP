function varargout=retop(varargin)
%  developpement en ondes planes dans un milieu homogene d'indice n 
% 		%%%%%%%%%%%%%%%%
% 		%    2 D       %
% 		%%%%%%%%%%%%%%%%
% [Ep,Em,angles]=retop(n,u,v,e,x,y,z,wx,wy,wz,k0,parm);
%
% 
% n:indice  du milieu homogene k0=2*pi/ld
% u,v: incidences ou on veut calculer le developpement en ondes planes
%            E = somme_sur_u_v de (tf(E) *exp(i (u*x+v*y)) )
% e(z,x,y,1:6) champ calcule par retchamp (6 composantes)
% wx,wy,wz:poids pour l'integration de Fourier
% l'un de ces vecteurs est vide :le developpement en ondes planes est fait dans la direction perpendiculaire
% parm=struct('apod',0,'delt',0,'uvmesh',0); (par defaut)
%     apod:pour l'integration de Fourier 
%     delt:pour la singularite de l'incidence normale 
% si uvmesh=1 u et v sont des vecteurs de meme dimension donnant les couples de valeurs de u et v
%  sinon u et v sont 2 vecteurs.Le calcul est fait sur ndgrid(u,v)
%    dans les tableaux suivants u ,v sont alors 'soudés' et ces tableaux perdent une dimension
%
% si wz=[] Ep(u,v,1:6,z),Em(u,v,1:6,z),sont les composantes de Fourier de e(z,x,y,1:6) par rapport à x,y fonctions de z (x->u y->v)
% si wx=[] Ep(u,v,1:6,x),Em(u,v,1:6,x),sont les composantes de Fourier de e(z,x,y,1:6) par rapport à y,z fonctions de x (y->u z->v)
% si wy=[] Ep(u,v,1:6,y),Em(u,v,1:6,y),sont les composantes de Fourier de e(z,x,y,1:6) par rapport à z,x fonctions de y (z->u x->v)
% ces composantes sont definies dans le meme repere que e (Ep(u,v,1,z) est Ex etc ...)
%*******************************************************************************
%  mise en forme des ondes planes en 2D (si angles est en sortie)
%
% si length(u)=nu,length(v)=nv,length(z)=nz:
% size(Ep)= size(Em)= [nu,nv,6,nz]
%          teta: [ nu nv ] angle de oz avec kp (0 a pi/2)
%         delta: [ nu nv ] angle de oy avec vp oriente par oz  (-pi a pi)
%          psip: [ nu nv nz] angle(oriente par kp) entre up et la direction principale de polarisation de Ep
%            kp: [ nu nv 3]  vecteur d'onde UNITAIRE dirige vers le haut
%            vp: [ nu nv 3]  unitaire perpendiculaire au plan d'incidence  (le triedre 0z, kp ,vp etant direct)  
%            up: [ nu nv 3]  up=vp /\ kp (le triedre up vp kp est direct)
%           EEp: [ nu nv 2 nz] composantes de Ep sur up et vp
%           HHp: [ nu nv 2 nz] composantes de Hp sur up et vp
%          EEEp: [ nu nv 2 nz] amplitude de Ep dans le repere deduit de up vp par rotation d'angle psip
%          HHHp: [ nu nv 2 nz] amplitude de Ep dans le repere deduit de up vp par rotation d'angle psip
%  angles.   Fp: [ nu nv nz]diagramme de diffraction vers le haut en energie (*1/r^2 ...)
%          psim: [ nu nv nz]:angle(oriente par km) entre um et la direction principale de polarisation de Em
%            km: [ nu nv 3] vecteur d'onde UNITAIRE dirige vers le bas
%            vm: [ nu nv 3] =vm unitaire perpendiculaire au plan d'incidence  (le triedre 0z, km ,vm est direct)  
%            um: [ nu nv 3] = vm /\ km (le triedre um vm km est direct) 
%           EEm: [ nu nv 2 nz] composantes de Em sur um et vm
%           HHm: [ nu nv 2 nz] composantes de Hm sur um et vm
%          EEEm: [ nu nv 2 nz] amplitude de Em dans le repere deduit de um vm par rotation d'angle psim
%          HHHm: [ nu nv 2 nz] amplitude de Hm dans le repere deduit de um vm par rotation d'angle psim
%            Fm: [ nu nv nz]diagramme de diffraction vers le bas en energie (*1/r^2 ...)
% pour l'incidence normale   (cas degenere) vp=vm=[-sin(delt),cos(delt)] 
% (l'angle delt peut etre introduit en parametre et vaut 0 par defaut)
% 
% 		%%%%%%%%%%%%%%%%
% 		%    1 D       %
% 		%%%%%%%%%%%%%%%%
% [Ep,Em,angles]=retop(n,u,e,x,y,wx,wy,pol,k0,parm);
%
% n:indice  du milieu homogene k0=2*pi/ld pol:0 E//  2 H//
% u: incidences où on veut calculer le developpement en ondes planes ( k0*n*sin(teta) )
%            E = somme_sur_u de (tf(E) *exp(i (u*x) )
% e(y,x,:) champ calcule par retchamp (3 ou 4 composantes)
% wx,wy:poids pour l'integration de Fourier
% l'un de ces vecteurs est vide :le developpement en ondes planes est fait dans la direction perpendiculaire
% parm=struct('apod',0); (par defaut)
%     apod:pour l'intégration de Fourier 
% 
% si wx=[] Ep(u,1:3,x),Em(u,1:3,x),sont les composantes de Fourier de e(y,x,1:3) par rapport à y fonctions de x
% si wy=[] Ep(u,1:3,y),Em(u,1:3,y),sont les composantes de Fourier de e(y,x,1:3) par rapport à x fonctions de y
%
%*******************************************************************************
%  mise en forme des ondes planes  en 1D (si angles est en sortie)
% si length(u)=nu,length(y)=ny:
% size(Ep)= size(Em)= [nu,3,nz]
%         teta: [ nu 1 ] angle de kp avec Oy orienté par Oz (-pi/2 a pi/2)
%           kp: [ nu 2 ] vecteur d'onde UNITAIRE dirige vers le haut
%           up: [ nu 2 ] perpendiculaire a kp  (Oz up kp direct) 
%          EEp: [ nu ny] composante de Ep sur Oz
%          HHp: [ nu ny] composante de Hp sur up
% angles.   Fp: [ nu ny] diagramme de diffraction vers le haut en energie (*1/r ...)
%           km: [ nu 2 ] vecteur d'onde UNITAIRE dirige vers le bas
%           um: [ nu 2 ] perpendiculaire a km  (Oz um km direct) 
%          EEm: [ nu ny] composante de Em sur Oz
%          HHm: [ nu ny] composante de Hm sur um
%           Fm: [ nu ny] diagramme de diffraction vers le bas en energie (*1/r ...)
%*******************************************************************************
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% correction de Haitao en 1D %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% on ajoute des plasmons sur les cotés.Ces plasmons doivent avoir été calculés auparavent par retpl (ou autrement)
% e doit etre remplacé par: {e,Plg,Pld} où : Pl={pl,x0,[nh,nb],y_dioptre}
% 	pl: amplitude du plasmon comme donné par retpl : donc normalisé par le flux du vecteur de poynting en x=0
%                pl est un scalaire (moyenne ou valeur en un point)
% 	x0: point à partir duquel on ajoute le plasmon
% 	nh: en haut
% 	nb: en bas
% 	y_dioptre: dioptre métal 
% Dans ce cas l'apodisation n'agit que sur le champ moins le plasmon
%
%
% 
% 		%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 		%      METHODE DE L'INTEGRALE ou DE LA BOITE        %
% 		%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  On integre Lorentz sur une boite à l'exterieur de laquelle le milieu est stratifié
%      un des champs est celui calculé ( qui doit satisfaire des COS)
%      l'autre est la réponse du milieu stratifié à une onde plane incidente
%       %%%%%%%%%%%%%%
% 		%     1 D    %
% 		%%%%%%%%%%%%%%
%   angles=retop(n,y,u,  e_x,x_x,y_x,w_x,   e_y,x_y,y_y,w_y, pol,k0,sens);
%
%  n: indices de haut en bas
%  y: cotes des dioptres de haut en bas ( méme origine que dans le calcul des champs)
% u: incidences où on veut calculer le developpement en ondes planes ( k0*n*sin(teta) )
% e_x,x_x,y_x,w_x champ calculé par retchamp sur 2 coupes en y ( y_x) w_x poids pour l'integration en x
% e_y,x_y,y_y,w_y champ calculé par retchamp sur 2 coupes en x ( x_y) w_y poids pour l'integration en y
% real(sens)  1 developpement au dessus   en y(1)    imag(sens)  0 les strates sont perpendiculaires à oy 
%            -1 developpement au dessous  en y(end)                 1 les strates sont perpendiculaires à ox 
%              y ne doit pas etre [] . Pour un milieu homogene n a 2 elements egaux
%  angles structure à champs: 'teta', 'k', 'u', 'EE', 'HH', 'F', 'Flux_poynting'
%
% Modes 
%-------
% angles=retop([n,neff],y,  e_x,x_x,y_x,w_x,   e_y,x_y,y_y,w_y, pol,k0,sens);
%        attention ne pas mettre pol après neff
% angles structure à champs de champs : amp_p, Fp=abs(amp_p)^2 (mode à droite)
%                                       amp_m, Fm=abs(amp_m)^2 (mode à gauche)
%                                       Flux_pornting ( sur la boite)
%    et aussi les champs 0D(non clonés) ayant servi à engendrer les modes 
%		
%       %%%%%%%%%%%%%%
% 		%     2 D    %
% 		%%%%%%%%%%%%%%
% angles=retop(n,z,u,v, e_xy,x_xy,y_xy,z_xy,wx_xy,wy_xy,   e_zx,x_zx,y_zx,z_zx,wz_zx,wx_zx, e_yz,x_yz,y_yz,z_yz,wy_yz,wz_yz,k0,sens,parm);  
% e_xy,x_xy,y_xy,z_xy,wx_xy,wy_xy, ...% champ sur 2 coupes perpendiculaires à oz avec les points et les poids
% e_zx,x_zx,y_zx,z_zx,wz_zx,wx_zx, ...% champ sur 2 coupes perpendiculaires à oy avec les points et les poids
% e_yz,x_yz,y_yz,z_yz,wy_yz,wz_yz,... % champ sur 2 coupes perpendiculaires à ox avec les points et les poids
% real(sens)= 1 developpement au dessus   en z(1)      imag(sens)=  0 les strates sont perpendiculaires à oz 
% real(sens)=-1 developpement au dessous  en z(end)    imag(sens)=  1 les strates sont perpendiculaires à ox 
%                                                      imag(sens)= -1 les strates sont perpendiculaires à oy 
%              z ne doit pas etre [] . Pour un milieu homogene n a 2 elements egaux
%  parm structure     defaultopt=struct('uvmesh',0,'clone',0,'test',0);
% si uvmesh=1 u et v sont des vecteurs de meme dimension donnant les couples de valeurs de u et v
%  angles structure à champs: 'teta', 'delta', 'psi', 'k', 'v', 'u' ,'EE' ,'HH' ,'EEE' ,'HHH' ,'F', 'Flux_poynting' 
%
% Modes cylindriques
%-----------------------
% angles=retop([n,neff,pol],z,LL,teta, e_xy,x_xy,y_xy,z_xy,wx_xy,wy_xy,   e_zx,x_zx,y_zx,z_zx,wz_zx,wx_zx, e_yz,x_yz,y_yz,z_yz,wy_yz,wz_yz,k0,sens); 
% imag(sens)  comme pour les ondes planes (seule la partie imaginaire sert pour dire la direction des strates)
% si LL=[] calcul direct pour les teta
% si LL est non vide  si LL=L1 LLL=L1,  si LL=[L1,L2] LLL=L1:L2,  si LL=[L1,L2,L3], LLL=L1:L2:L3, sinon LL=retelimine de LL
% si LL est un vecteur colonne LLL=LL.'
% attention avec , si LL=-1:1 on obtient LLL=-1:0:1 =[] donc il vaut mieux faire LL=(-1:1).' 
%      calcul des coefficients de Fourier LLL, puis synthese de Fourier en teta
% angles structure à champs de champs f, amp, F f coefficients de Fourier, amp amplitude(teta) F intensite(teta), LLL
%    et aussi les champs 0D(non clonés)ayant servi à engendrer les modes. 
% si [teta,wteta]=retgauss(0,2*pi,...), l'energie sur le mode est sum(angles.F.*wteta(:)) ou 2*pi*sum(abs(f).^2)
%
% Modes de Choon
%-----------------------
% angles=retop( e_mode, x_mode,z_mode,neff, e_xy,x_xy,y_xy,z_xy,wx_xy,wy_xy,   e_zx,x_zx,y_zx,z_zx,wz_zx,wx_zx, e_yz,x_yz,y_yz,z_yz,wy_yz,wz_yz,k0,sens,parm); 
%     ------------------------------------------
%     | sens | propagation | x_mode | z_mode   |
%     ------------------------------------------
%     |  1   |   Oy       |  x     |  z       |
%     |  -1  |   Ox       |  y     |  z       |
%     | 1+i  |   Oz       |  y     |  x       |
%     | -1+i |   Oy       |  z     |  x       |
%     | 1-i  |   Ox       |  z     |  y       |
%     |-1-i  |   Oz       |  x     |  y       |
%     -----------------------------------------
% e_mode, x_mode,z_mode,neff calculés par retmode neff=gama/k0
%  parm structure à champs facultatif ( si parm.test=1, impression des champs  sur la boite )
% 
% 
%
%       %%%%%%%%%%%%%%%%%%%%%%%%%%%
% 		%    cylindres Popov      %
% 		%%%%%%%%%%%%%%%%%%%%%%%%%%%
% angles=retop(n,z,u,v,  e_x,x_x,z_x,wx_x,    e_z,x_z,z_z,wz_z,   k0,sens,L,sym,parm);    
% e_x,x_x,z_x,wx_x      champ sur 2 coupes perpendiculaires à oz x_x variant de 0 à r avec les points et les poids
% e_z,x_z,z_z,wz_z      champ pour x_z=r 
% TRES IMPORTANT:
%  Le calcul des champs doit etre IMPERATIVEMENT effectué avec l'option sym =0 meme si on calcule le developpement avec sym=1 ou -1
% ( recalculer init ou faire init{end}.sym=0 avant retchamp )
%
%  sens  1 developpement au dessus   en z(1)
%       -1 developpement au dessous  en z(end)
%  L sym parametres de popov
%
%  parm,u,v,angles:  comme en 2D
%
% Modes cylindriques
%-----------------------
% angles=retop([n,neff,pol],z,L,teta,   e_x,x_x,z_x,wx_x,    e_z,x_z,z_z,wz_z, k0,sym);
% angles structure à champs de champs f, amp, F f coefficient de Fourier sur L, amp amplitude(teta) tenant compte du parametre sym 
% F intensite(teta)  tenant compte du parametre sym et aussi les champs 0D(non clonés)ayant servi à engendrer les modes   
%--------------------------------------------------------------------------
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     %  AUTOMATISATION DE LA METHODE DE LA BOITE  %
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% INITIALISATION
%
% init=retop(cal_champ,k0,n_strates,z_strates,Centre,L,N_gauss,   parm); en 2D
% init=retop(cal_champ,pol,k0,n_strates,z_strates,Centre,L,N_gauss,   parm); en 1D
% init=retop(cal_champ,k0,n_strates,z_strates,L_popov,z_Centre,HR,N_gauss, parm); cylindres Popov
%
%   
% parm=struct('orientation',0,'disc',{{[],[],[]}},'option_cal_champ',2,'option_SI',0); par defaut
%     si option_cal_champ=0:
%     cal_champ doit etre une fonction qui calcule le champ non cloné au format retchamp 
%     e =calchamp(x(:),y(:),z(:))    size e=[length(z(:)) , length(x(:)) , length(y(:)) , 6] en 2D
%     e =calchamp(x(:),y(:))    size e=[length(y(:)) , length(x(:)) , 3] en 1D
%     si option_cal_champ=1 ou 2:
%     cal_champ doit etre une fonction qui calcule le champ non cloné au format suivant
%     e =calchamp(x(:),y(:),z(:))    size e=[length(x(:))=length(y(:))=length(z(:)) , 6] en 2D
%     e =calchamp(x(:),y(:))    size e=[length(x(:))=length(y(:)) , 3] en 1D
%     si option_cal_champ=2 la fonction cal_champ n'est appeleé qu'une fois 
% option_SI<0 exp(-iwt) option_SI>0 exp(iwt)
% Scale=abs(option_SI) unite de longueur en metres (ld x y z k0 u v)
% si option_SI~=0 E et H sont normalisés avec le facteur Z0 
%   E*sqrt(Z0) H/sqrt(Z0)
% option_SI=1 comsol
% option_SI=-1.e6 conventions reticolo
% option_SI=0 conventions reticolo
%  
% k0=2*pi/ld;
% n_strates,z_strates,  length(n_strates) = 1+ length(z_strates) >2  de haut en bas
% z_strates ne peut par etre []
% L'ORIGINE du DOP est z_strates(1) si sens=1 (en haut)
%                      z_strates(end) si sens=-1 (en bas)
% Centre:Centre de la boite (longueur 3 en 2D, 2 en 1D)
% L  largeurs de la boite
% N_gauss pour points de gauss tableau de size [3,2]
% par defaut:N_gauss=[[10,10];[10,10];[10,10]])
% premiere colonne: degre, seconde repetition
% La seconde colonne est facultative
%
%   orientation, disc facultatifs
% en 2D:
%     orientation:  0 (par defaut) les strates sont perpendiculaires à oz 
%                   1 les strates sont perpendiculaires à ox 
%                  -1 les strates sont perpendiculaires à oy 
% en 1D:
%     orientation:  0 (par defaut) les strates sont perpendiculaires à oy 
%                   1 les strates sont perpendiculaires à ox 
% disc cell array des discontinuites pour gauss  
%
% CALCUL DU DOP
% en 2D:
% angles=retop(init,u,v,sens,  parm );
% u v: K//
% sens=1 au dessus, -1 au dessous
% parm=struct('uvmesh',1,'test',0); par defaut facultatif  test=1: tracé de la boite (recommandé)
% Attention l'option uvmesh est l'inverse de la version classique: on parametrise plus souvent en teta et phi que en u et v
% en 1D: angles=retop(init,u,sens,parm);parm=struct('test',0); par defaut
%
% Popov: angles=retop(init,u,v,sens,parm);
% parm=struct('uvmesh',1,'sym',0,'test',0); par defaut
% Pour les modes: en 2D
% angles=retop(init,neff,pol,L,teta, parm); pol=0 ou 2 polarisation du mode
% L=-5:5 ...teta est obligatoire et peut etre []  
% en 1D:angles=retop(init,neff, parm);parm=struct('test',0) par defaut. test=1 tracé de la boite
% popov: angles=retop(init,neff,pol,parm);
% parm=struct('sym',0,'test',0,'teta',[]); par defaut
% si [teta,wteta]=retgauss(0,2*pi,...), l'energie sur le mode est sum(angles.F.*wteta(:)) ou 2*pi*sum(abs(f).^2)
%--------------------------------------------------------------------------
%
% % EXEMPLE 2 D
% n=1.5;k0=2*pi/1.2;% Fourier en x z
% y=[-2,2]/(k0*n);wy=[];[z,wz]=retgauss(-100/(k0*n),100/(k0*n),10,25,[-6,6]*max(abs(y)));[x,wx]=retgauss(-100/(k0*n),100/(k0*n),10,25,[-6,6]*max(abs(y)));
% pols=rand(1,6)+i*rand(1,6);e=retpoint(k0*n^2,k0,pols,[0,0,0],x,y,z);% source ponctuelle
% [teta,wt]=retgauss(0,pi/2,25);[delta,wd]=retgauss(0,2*pi,-50);nt=length(teta);nd=length(delta);
% [Teta,Delta]=ndgrid(teta,delta);[wt,wd]=ndgrid(wt,wd);w=wt(:).*wd(:).*sin(Teta(:));
% u=k0*n*sin(Teta).*cos(Delta);v=k0*n*sin(Teta).*sin(Delta);
% [ep,em,angles]=retop(n,u(:),v(:),e,x,y,z,wx,wy,wz,k0,struct('apod',7.25,'uvmesh',1));
%  energie_vers_le_haut=w.'*angles.Fp(:,2)
%  energie_vers_le_bas=w.'*angles.Fm(:,1)
%
% %  EXEMPLE 1 D
% %Fourier en x
% n=1.5;k0=2*pi/1.2;pol=2;
% y=[-2,2]/(k0*n);wy=[];[x,wx]=retgauss(-300/(k0*n),300/(k0*n),20,20,[-6,6]*max(abs(y)));
% pols=randn(1,3)+i*randn(1,3);
% e=retpoint(k0*n^2,k0,pols,[0,0],x,y,pol);
% [teta,wt]=retgauss(-pi/2,pi/2,20,20);u=k0*n*sin(teta);
% [ep,em,angles]=retop(n,u,e,x,y,wx,wy,pol,k0,struct('apod',7.25));%figure;plot(u,squeeze(abs(ep)).^2,u,squeeze(abs(em)).^2,'.')
% Ap=angles.Fp.'*wt.';
% Am=angles.Fm.'*wt.';
% pt=retpoynting(e,[0,1]);Pt=diag([-1,1])*pt*wx.';
% retcompare(Pt,[Am(1),Ap(2)])
%
%
% See also:RETB,RETCHAMP,RETAPOD,RETGAUSS,RETPL,RETHELP_POPOV,RET_NORMALISE_DOP
%


%%%%%%%%%%%%%%%%%%%%
%   AIGUILLAGE     %
%%%%%%%%%%%%%%%%%%%%



if isa(varargin{1},'function_handle');% methode de la boite construction automatique et calcul du champ
switch  length(varargin{6})   
case 3;[varargout{1:nargout}]=boite_automatique_init_2D(varargin{:});
case 2;[varargout{1:nargout}]=boite_automatique_init_1D(varargin{:});
case 1;[varargout{1:nargout}]=boite_automatique_init_popov(varargin{:});
end;
return;end;

if iscell(varargin{1});
if nargin>8;[varargout{1:nargout}]=ondes_planes_popov(varargin{:});% cylindres Popov
else;
switch  length(varargin{1}{end})   
case 8;[varargout{1:nargout}]=boite_automatique_1D(varargin{:});   
case 18;[varargout{1:nargout}]=boite_automatique_2D(varargin{:});
case 9;[varargout{1:nargout}]=boite_automatique_popov(varargin{:});
end;
end
return;end;

if nargin>23;[varargout{1:nargout}]=caldop_2D(varargin{:});return;end;    % 2 D integrale (version 2010)
if nargin>16;[varargout{1:nargout}]=caldop_popov(varargin{:});return;end; % popov integrale (version 2010)
if nargin>12;
if length(varargin{1})>length(varargin{2})+2;[varargout{1:nargout}]=caldop_modes_popov(varargin{:});% modes circulaires symetrie popov
else;[varargout{1:nargout}]=caldop_1D(varargin{:}); % 1 D integrale (version 2010)
end;
return;end;
if nargin>10;[varargout{1:nargout}]=ondes_planes_2D(varargin{:});return;end;% 2 D
if nargin>5;[varargout{1:nargout}]=ondes_planes_1D(varargin{:});return;end;% 1 D
[varargout{1:nargout}]=ondes_planes_1D_ancien(varargin{:}); % compatibilitee avec une ancienne version 1 D

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function varargout=boite_automatique_2D(varargin);
if isstruct(varargin{end}) | isempty(varargin{end});prv=1;else prv=0;end
if nargin==4+prv & ~isempty(varargin{end}) ;[varargout{1:nargout}]=boite_automatique_op_2D(varargin{:});
else;[varargout{1:nargout}]=boite_automatique_modes_2D(varargin{:});
end;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function e=comsol2ret(e,option_SI,pol);% e comsol -> e ret
if option_SI>0;e=conj(e);end;
ratioH=19.409541814833513;%ratioH=sqrt(retconstantes('Z0'));
ratioE=51.5210513231055e-3;%ratioE=1./sqrt(retconstantes('Z0'));
switch(pol);
case 1;e(:,1:3)=e(:,1:3)*ratioE;e(:,4:6)=e(:,4:6)*ratioH;% 3D 
case 0;e(:,1)=e(:,1)*ratioE;e(:,2:3)=e(:,2:3)*ratioH;% TE 
case 2;e(:,1)=e(:,1)*ratioH;e(:,2:3)=-e(:,2:3)*ratioE;% TM 
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function init_dop=boite_automatique_init_2D(cal_champ,k0,n_strates,z_strates,Centre,L,N_gauss,parm);
defaultopt=struct('orientation',0,'disc',{{[],[],[]}},'option_cal_champ',2,'option_SI',0,'option_i',nan);
% definition de la boite et calcul du champ
if nargin<8;parm=[];end;
if nargin<7;N_gauss=[[10,10];[10,10];[10,10]];end;if numel(N_gauss)==3;N_gauss=[[10;10;10],ceil(N_gauss(:)/10)];end;
orientation=retoptimget(parm,'orientation',defaultopt,'fast');
disc=retoptimget(parm,'disc',defaultopt,'fast');
option_cal_champ=retoptimget(parm,'option_cal_champ',defaultopt,'fast');
option_SI=retoptimget(parm,'option_SI',defaultopt,'fast');
option_i=retoptimget(parm,'option_i',defaultopt,'fast');
if ~isnan(option_i);option_SI=option_i;end;
if isempty(z_strates);z_strates=Centre(3);n_strates=[n_strates,n_strates];end;
if option_SI~=0;% passage m-> microns meme pour k0
Scale=1.e-6/abs(option_SI);
z_strates=z_strates/Scale;Centre=Centre/Scale;L=L/Scale;
k0=k0*Scale;
cal_champ=@(x,y,z) comsol2ret(cal_champ(x*Scale,y*Scale,z*Scale),option_SI,1); 
end;
x_yz=[Centre(1)-L(1)/2,Centre(1)+L(1)/2];
y_zx=[Centre(2)-L(2)/2,Centre(2)+L(2)/2];
z_xy=[Centre(3)-L(3)/2,Centre(3)+L(3)/2];
switch orientation;% tres important
case 1;disc{1}=[disc{1}(:).',z_strates(:).'];
case -1;disc{2}=[disc{2}(:).',z_strates(:).'];
case 0;disc{3}=[disc{3}(:).',z_strates(:).'];
end;

[x_zx,wx_zx]=retgauss(x_yz(1),x_yz(2),N_gauss(1,1),N_gauss(1,2),disc{1});
[y_yz,wy_yz]=retgauss(y_zx(1),y_zx(2),N_gauss(2,1),N_gauss(2,2),disc{2});
[z_zx,wz_zx]=retgauss(z_xy(1),z_xy(2),N_gauss(3,1),N_gauss(3,2),disc{3});
x_xy=x_zx;wx_xy=wx_zx;y_xy=y_yz;wy_xy=wy_yz;z_yz=z_zx;wz_yz=wz_zx;

if option_cal_champ~=0;
[Z_xy,X_xy,Y_xy]=ndgrid(z_xy,x_xy,y_xy);
[Z_zx,X_zx,Y_zx]=ndgrid(z_zx,x_zx,y_zx);
[Z_yz,X_yz,Y_yz]=ndgrid(z_yz,x_yz,y_yz);
n_xy=length(Z_xy(:));n_zx=length(Z_zx(:));n_yz=length(Z_yz(:));
end;
switch option_cal_champ;
case 0;% 3 appels size(e)=[nz, nx, ny,6] reticolo retchamp
e_xy=cal_champ(x_xy,y_xy,z_xy);
e_zx=cal_champ(x_zx,y_zx,z_zx);
e_yz=cal_champ(x_yz,y_yz,z_yz);
case 2;% un seul appel size(e)=[nbpoints,6]
Z=[Z_xy(:);Z_zx(:);Z_yz(:)];X=[X_xy(:);X_zx(:);X_yz(:)];Y=[Y_xy(:);Y_zx(:);Y_yz(:)];
e=cal_champ(X,Y,Z);
e_xy=reshape(e(1:n_xy,:),length(z_xy),length(x_xy),length(y_xy),6);
e_zx=reshape(e(1+n_xy:n_zx+n_xy,:),length(z_zx),length(x_zx),length(y_zx),6);
e_yz=reshape(e(1+n_zx+n_xy:n_yz+n_zx+n_xy,:),length(z_yz),length(x_yz),length(y_yz),6);
case 1;% 3 appels size(e)=[nbpoints,6]
e_xy=reshape(cal_champ(X_xy,Y_xy,Z_xy),length(z_xy),length(x_xy),length(y_xy),6);
e_zx=reshape(cal_champ(X_zx,Y_zx,Z_zx),length(z_zx),length(x_zx),length(y_zx),6);
e_yz=reshape(cal_champ(X_yz,Y_yz,Z_yz),length(z_yz),length(x_yz),length(y_yz),6);
end;
champ={ e_xy,x_xy,y_xy,z_xy,wx_xy,wy_xy,   e_zx,x_zx,y_zx,z_zx,wz_zx,wx_zx, e_yz,x_yz,y_yz,z_yz,wy_yz,wz_yz};
init_dop={n_strates,z_strates,k0,orientation+1i*option_SI,champ};

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function angles=boite_automatique_op_2D(init,u,v,sens,parm);
defaultopt=struct('option_parfor',0,'uvmesh',1,'test',0);
if nargin<5;parm=struct('');end;if isempty(parm);parm=struct('');end;
uvmesh=retoptimget(parm,'uvmesh',defaultopt,'fast');
option_parfor=retoptimget(parm,'option_parfor',defaultopt,'fast');
test=retoptimget(parm,'test',defaultopt,'fast');
[n_strates,z_strates,k0,orientation,champ]=deal(init{:});
option_SI=imag(orientation);orientation=real(orientation);
if option_SI~=0;Scale=1.e-6/abs(option_SI);u=u*Scale;v=v*Scale;else;Scale=1;end;
angles=retop(n_strates,z_strates,u(:),v(:),champ{:},k0,sens+1i*orientation,struct('option_parfor',option_parfor,'uvmesh',uvmesh,'test',test,'Scale',Scale));

if option_SI~=0;% mise en forme et conversion en unite SI
angles.EH=angles.E;    
angles=rmfield(angles,{'e','E','teta','Flux_poynting'});

ratioH=19.409541814833513; % ratioH=sqrt(retconstantes('Z0'));
ratioE=51.5210513231055e-3;% ratioE=1./sqrt(retconstantes('Z0'));
if option_SI>0;
angles.EH=conj(angles.EH);angles.EE=conj(angles.EE);%angles.EEE=conj(angles.EEE);
angles.HH=conj(angles.HH);angles.HHH=conj(angles.HHH);
end;
angles.EH(:,1:3)=angles.EH(:,1:3)/ratioE;
angles.EH(:,4:6)=angles.EH(:,4:6)/ratioH;
angles.EE=angles.EE/ratioE;angles.EEE=angles.EEE/ratioE;
angles.HH=angles.HH/ratioH;angles.HHH=angles.HHH/ratioH;
angles.origine=angles.origine*Scale;

angles.F=angles.F*Scale^2;   % modif 11 2016
angles.EE=angles.EE*Scale^2; % modif 11 2016
angles.HH=angles.HH*Scale^2; % modif 11 2016
angles.EH=angles.EH*Scale; % modif 11 2016

angles=orderfields(angles,[1,2,3,4,5,6,13,7,8,9,10,11,12]);
angles=rmfield(angles,{'EEE','HHH','psi'});

end;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function angles=boite_automatique_modes_2D(init,neff,pol,LL,teta,parm);
if nargin<6;parm=struct('');end;if isempty(parm);parm=struct('');end;
[n_strates,z_strates,k0,orientation,champ]=deal(init{:});
option_SI=imag(orientation);orientation=real(orientation);
if option_SI~=0;Scale=1.e-6/abs(option_SI);else;Scale=1;end;if ~isempty(parm);parm.Scale=Scale;end;
n_strates=[n_strates,neff,pol];
angles=retop(n_strates,z_strates,LL(:),teta,champ{:},k0,1i*orientation,parm);
if option_SI~=0;% mise en forme et conversion en unite SI
angles.F=angles.F*Scale^2;angles.amp=angles.amp*Scale; % modif 11 2016
angles=rmfield(angles,{'Flux_poynting','exy','ezx','eyz'});
end;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function init_dop=boite_automatique_init_1D(cal_champ,pol,k0,n_strates,z_strates,Centre,L,N_gauss,parm);
defaultopt=struct('orientation',0,'disc',{{[],[]}},'option_cal_champ',2,'option_SI',0,'option_i',nan);
% definition de la boite et calcul du champ
if nargin<9;parm=struct('');end;if isempty(parm);parm=struct('');end;
if nargin<8;N_gauss=[[10,10];[10,10]];end;if numel(N_gauss)==2;N_gauss=[[10;10],ceil(N_gauss(:)/10)];end;
orientation=retoptimget(parm,'orientation',defaultopt,'fast');
disc=retoptimget(parm,'disc',defaultopt,'fast');
option_cal_champ=retoptimget(parm,'option_cal_champ',defaultopt,'fast');
option_SI=retoptimget(parm,'option_SI',defaultopt,'fast');
option_i=retoptimget(parm,'option_i',defaultopt,'fast');
if ~isnan(option_i);option_SI=option_i;end;
if isempty(z_strates);z_strates=Centre(2);n_strates=[n_strates,n_strates];end;
if option_SI~=0;% passage m-> microns meme pour  k0
Scale=1.e-6/abs(option_SI);
z_strates=z_strates/Scale;Centre=Centre/Scale;L=L/Scale;
k0=k0*Scale;
cal_champ=@(x,y) comsol2ret(cal_champ(x*Scale,y*Scale),option_SI,pol); 
end;

switch orientation;
case 1;disc{1}=[disc{1}(:).',z_strates(:).'];
case 0;disc{2}=[disc{2}(:).',z_strates(:).'];
end;
x_y=[Centre(1)-L(1)/2,Centre(1)+L(1)/2];
y_x=[Centre(2)-L(2)/2,Centre(2)+L(2)/2];
[x_x,w_x]=retgauss(x_y(1),x_y(2),N_gauss(1,1),N_gauss(1,2),disc{1});
[y_y,w_y]=retgauss(y_x(1),y_x(2),N_gauss(2,1),N_gauss(2,2),disc{2});

if option_cal_champ~=0;
[Y_x,X_x]=ndgrid(y_x,x_x);
[Y_y,X_y]=ndgrid(y_y,x_y);
n_x=length(Y_x(:));n_y=length(Y_y(:));
end;

switch option_cal_champ;
case 0;% 2 appels size(e)=[nz, nx, ny,6] reticolo retchamp
e_x=cal_champ(x_x,y_x);
e_y=cal_champ(x_y,y_y);
case 1;% 2 appels size(e)=[nbpoints,6]
e_x=reshape(cal_champ(X_x,Y_x),length(y_x),length(x_x),3);
e_y=reshape(cal_champ(X_y,Y_y),length(y_y),length(x_y),3);
case 2;% un seul appel size(e)=[nbpoints,6]
Y=[Y_x(:);Y_y(:)];X=[X_x(:);X_y(:)];
e=cal_champ(X,Y);
e_x=reshape(e(1:n_x,:),length(y_x),length(x_x),3);
e_y=reshape(e(1+n_x:n_x+n_y,:),length(y_y),length(x_y),3);
end;


champ={ e_x,x_x,y_x,w_x,   e_y,x_y,y_y,w_y};
init_dop={n_strates,z_strates,pol,k0,orientation+1i*option_SI,champ};

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function varargout=boite_automatique_1D(varargin);
if isstruct(varargin{end}) | isempty(varargin{end});prv=1;else prv=0;end
if nargin==(3+prv);[varargout{1:nargout}]=boite_automatique_op_1D(varargin{:});
else;[varargout{1:nargout}]=boite_automatique_modes_1D(varargin{:});end;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function angles=boite_automatique_op_1D(init,u,sens,parm);
if nargin<4;parm=struct('');end;if isempty(parm);parm=struct('');end;
[n_strates,z_strates,pol,k0,orientation,champ]=deal(init{:});
option_SI=imag(orientation);orientation=real(orientation);
if option_SI~=0;Scale=1.e-6/abs(option_SI);u=u*Scale;else;Scale=1;end;if ~isempty(parm);parm.Scale=Scale;end;
angles=retop(n_strates,z_strates,u(:),champ{:},pol,k0,sens+1i*orientation,parm);
if option_SI~=0;% mise en forme et conversion en unite SI
angles.EH=angles.e;
angles=rmfield(angles,{'e','Flux_poynting'});

ratioH=19.409541814833513; % ratioH=sqrt(retconstantes('Z0'));
ratioE=51.5210513231055e-3;% ratioE=1./sqrt(retconstantes('Z0'));
if option_SI>0;
angles.EH=conj(angles.EH);angles.EE=conj(angles.EE);angles.HH=conj(angles.HH);
end;

switch(pol);
case 0;
angles.EH(:,1)=angles.EH(:,1)/ratioE;
angles.EH(:,2:3)=angles.EH(:,2:3)/ratioH;
[angles.EE,angles.HH]=deal(angles.EE/ratioE,angles.HH/ratioH);
case 2;
angles.EH(:,1)=angles.EH(:,1)/ratioH;
angles.EH(:,2:3)=-angles.EH(:,2:3)/ratioE;
[angles.HH,angles.EE]=deal(angles.EE/ratioH,-angles.HH/ratioE);
end;

angles.origine=angles.origine*Scale;
angles.F=angles.F*Scale;% modif 11 2016
angles.EH=angles.EH*sqrt(Scale);% modif 11 2016
angles.EE=angles.EE*Scale;% modif 11 2016
angles.HH=angles.HH*Scale;% modif 11 2016
angles=orderfields(angles,[1,2,3,8,4,5,6,7]);
end;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function angles=boite_automatique_modes_1D(init,neff,parm);
if nargin<3;parm=struct('');end;if isempty(parm);parm=struct('');end;
[n_strates,z_strates,pol,k0,orientation,champ]=deal(init{:});
option_SI=imag(orientation);orientation=real(orientation);
n_strates=[n_strates,neff];
if option_SI~=0;Scale=1.e-6/abs(option_SI);else;Scale=1;end;if ~isempty(parm);parm.Scale=Scale;end;
angles=retop(n_strates,z_strates,champ{:},pol,k0,1i*orientation,parm);
if option_SI~=0;% mise en forme et conversion en unite SI
angles=rmfield(angles,{'ex','ey','Flux_poynting'});
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function init_dop=boite_automatique_init_popov(cal_champ,k0,n_strates,z_strates,L,z_Centre,HR,N_gauss,parm); %cylindres Popov
% angles=retop(n,z,u,v,  e_x,x_x,z_x,wx_x,    e_z,x_z,z_z,wz_z,   k0,sens,L,sym,parm);    
% e_x,x_x,z_x,wx_x      champ sur 2 coupes perpendiculaires à oz x_x variant de 0 à r avec les points et les poids
% e_z,x_z,z_z,wz_z      champ pour x_z=r 

defaultopt=struct('disc',{{[],[]}},'option_cal_champ',2);
% definition de la boite et calcul du champ
if nargin<9;parm=struct('');end;if isempty(parm);parm=struct('');end;
if nargin<8;N_gauss=[[10,10];[10,10]];end;
disc=retoptimget(parm,'disc',defaultopt,'fast');
option_cal_champ=retoptimget(parm,'option_cal_champ',defaultopt,'fast');

disc{1}=[disc{1}(:).',z_strates(:).'];

z_x=[z_Centre-HR(1)/2,z_Centre+HR(1)/2];
x_z=HR(2);
[z_z,wz_z]=retgauss(z_x(1),z_x(2),N_gauss(1,1),N_gauss(1,2),disc{1});
[x_x,wx_x]=retgauss(0,x_z,N_gauss(2,1),N_gauss(2,2),disc{2});

if option_cal_champ~=0;
[Z_x,X_x,Y_x]=ndgrid(z_x,x_x,0);
[Z_z,X_z,Y_z]=ndgrid(z_z,x_z,0);
n_x=length(Z_x(:));n_z=length(Z_z(:));
end;
switch option_cal_champ;
case 0;% 2 appels size(e)=[nz, nx, ny,6] reticolo retchamp
e_x=cal_champ(x_x,0,z_x);
e_z=cal_champ(x_z,0,z_z);
case 2;% un seul appel size(e)=[nbpoints,6]
Z=[Z_x(:);Z_z(:)];X=[X_x(:);X_z(:)];Y=[Y_x(:);Y_z(:)];
e=cal_champ(X,Y,Z);
e_x=reshape(e(1:n_x,:),length(z_x),length(x_x),1,6);
e_z=reshape(e(1+n_x:n_z+n_x,:),length(z_z),length(x_z),1,6);
case 1;% 2 appels size(e)=[nbpoints,6]
e_x=reshape(cal_champ(X_x,Y_x,Z_x),length(z_x),length(x_x),1,6);
e_z=reshape(cal_champ(X_z,Y_z,Z_z),length(z_z),length(x_z),1,6);
end;

champ={ e_x,x_x,z_x,wx_x,    e_z,x_z,z_z,wz_z  , k0};
init_dop={n_strates,z_strates,L,champ};
% retop(n,z,u,v,  e_x,x_x,z_x,wx_x,    e_z,x_z,z_z,wz_z,   k0,sens,L,sym,parm); 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function varargout=boite_automatique_popov(varargin);
if isstruct(varargin{end}) | isempty(varargin{end});nargin_nv=nargin;else nargin_nv=nargin+1;end;
if nargin_nv==5;[varargout{1:nargout}]=boite_automatique_op_popov(varargin{:});
else;[varargout{1:nargout}]=boite_automatique_modes_popov(varargin{:});
end;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function angles=boite_automatique_op_popov(init,u,v,sens,parm);
defaultopt=struct('uvmesh',1,'sym',0,'test',0);
if nargin<5;parm=struct('');end;if isempty(parm);parm=struct('');end;
uvmesh=retoptimget(parm,'uvmesh',defaultopt,'fast');
sym=retoptimget(parm,'sym',defaultopt,'fast');
test=retoptimget(parm,'test',defaultopt,'fast');
[n_strates,z_strates,L,champ]=deal(init{:});
angles=retop(n_strates,z_strates,u(:),v(:),champ{:},sens,L,sym,struct('uvmesh',uvmesh,'test',test));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function angles=boite_automatique_modes_popov(init,neff,pol,parm);
defaultopt=struct('sym',0,'test',0,'teta',[]);
if nargin<4;parm=struct('');end;if isempty(parm);parm=struct('');end;
sym=retoptimget(parm,'sym',defaultopt,'fast');
test=retoptimget(parm,'test',defaultopt,'fast');
teta=retoptimget(parm,'teta',defaultopt,'fast');
[n_strates,z_strates,L,champ]=deal(init{:});
n_strates=[n_strates,neff,pol];
angles=retop(n_strates,z_strates,L,teta,champ{:},sym,struct('test',test));


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function varargout=ondes_planes_1D(n,u,e,x,y,wx,wy,pol,k0,parm);
if nargin<10;parm=[];end;
if nargin<9;k0=1;end;
if nargin<8;pol=0;end;
if nargin<7;wy=[];end;
ep=retep(n,pol,k0);

defaultopt=struct('apod',0);
apod=retoptimget(parm,'apod',defaultopt,'fast');
cal_angles=nargout>2;
if iscell(e);e{1}=e{1}(:,:,1:3);else;e=e(:,:,1:3);end;% cas de la correction de Haitao


if isempty(wx);%[ x y ] -->[ y x]
if iscell(e);e{1}=permute(e{1},[2,1,3]);e{1}=e{1}(:,:,[1,3,2]);e{1}(:,:,1)=-e{1}(:,:,1); % cas de la correction de Haitao
else;e=permute(e,[2,1,3]);e=e(:,:,[1,3,2]);e(:,:,1)=-e(:,:,1);end;
[varargout{1:nargout}]=ondes_planes_1Dy(ep,u,e,y,x,wy,apod,cal_angles,k0);
for ii=1:min(2,nargout);varargout{ii}(:,1,:)=-varargout{ii}(:,1,:);varargout{ii}=varargout{ii}(:,[1,3,2],:);end;
return;end;
if isempty(wy);
[varargout{1:nargout}]=ondes_planes_1Dy(ep,u,e,x,y,wx,apod,cal_angles,k0);
return;end;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Ep,Em,angles]=ondes_planes_1Dy(ep,u,e,x,y,wx,apod,cal_angles,k0);
if iscell(e);% cas de la correction de Haitao
Plg=e{2};Pld=e{3};Haitao=1;
e=e{1};	
else;Haitao=0;end;
e(:,:,2:3)=i*e(:,:,2:3);% clonnage
nx=length(x);ny=length(y);nu=length(u);ne=3;
Ep=zeros(ny,nu,ne);Em=zeros(ny,nu,ne);

if Haitao;% on ajoute la Tf des plasmons sur les cotès
if apod~=0;epl=zeros(size(e));x1=min(x);x2=max(x);[fg,fd]=retfind(x<(x1+x2)/2);end;
if ~isempty(Pld);% droite
[pl_g,pl_d,eg,ed]=retpl(Pld{3},Pld{4},zeros(length(y),1,3),0,y,[],1,k0);	
ed(:,:,2:3)=i*ed(:,:,2:3);% clonage
if length(Pld{3})==2;beta_pl=retplasmon(Pld{3}(1),Pld{3}(2),2*pi/k0);beta_pl=beta_pl.constante_de_propagation;else;beta_pl=retmode(Pld{3}(end),Pld{3}(1:end-2),-diff(Pld{4})*k0,Pld{3}(end-1));end;
prv=i*(k0*beta_pl-u(:).');
for ie=1:2;
Ep(:,:,ie)=Ep(:,:,ie)-(Pld{1}/(2*pi))*(ed(:,1,ie)*(exp(prv*Pld{2})./prv));
if apod~=0;epl(:,fd,ie)=epl(:,fd,ie)+Pld{1}*ed(:,1,ie)*exp(i*k0*beta_pl*x(fd));end;
end;
end;	
if ~isempty(Plg); % gauche
[pl_g,pl_d,eg,ed]=retpl(Plg{3},Plg{4},zeros(length(y),1,3),0,y,[],1,k0);	
eg(:,:,2:3)=i*eg(:,:,2:3); % clonage
if length(Plg{3})==2;beta_pl=retplasmon(Plg{3}(1),Plg{3}(2),2*pi/k0);beta_pl=beta_pl.constante_de_propagation;else;beta_pl=retmode(Plg{3}(end),Plg{3}(1:end-2),-diff(Plg{4})*k0,Plg{3}(end-1));end;
prv=-i*(k0*beta_pl+u(:).');
for ie=1:2;
Ep(:,:,ie)=Ep(:,:,ie)+(Plg{1}/(2*pi))*(eg(:,1,ie)*(exp(prv*Plg{2})./prv));
if apod~=0;epl(:,fg,ie)=epl(:,fg,ie)+Plg{1}*eg(:,1,ie)*exp(-i*k0*beta_pl*x(fg));end;
end;
end;
end;% Haitao

%ax=retdiag(wx.*retchamp([apod,nx]))*exp(-i*x(:)*u(:).')/(2*pi);
if apod~=0;apodise=interp1(linspace(min(x),max(x),100),retchamp([apod,100]),x);
if Haitao;for ie=1:2;e(:,:,ie)=epl(:,:,ie)+(e(:,:,ie)-epl(:,:,ie))*retdiag(apodise);end;% on n'apodise que la partie sans plasmons
else;for ie=1:2;e(:,:,ie)=e(:,:,ie)*retdiag(apodise);end;end;
end;

ax=retdiag(wx)*exp(-i*x(:)*u(:).')/(2*pi);
for ie=1:2;Ep(:,:,ie)=Ep(:,:,ie)+e(:,:,ie)*ax;end;

khi=retsqrt(ep(1)*ep(2)-u.^2,-1);khi=repmat(khi,ny,1);
[Em(:,:,1),Em(:,:,2)]=deal(.5*(Ep(:,:,1)+(i/ep(3))*Ep(:,:,2)./khi),.5*(-(i*ep(3))*Ep(:,:,1).*khi+Ep(:,:,2)));
Em(:,:,3)=-(i/ep(2))*repmat(u,ny,1).*Em(:,:,1);
[Ep(:,:,1),Ep(:,:,2)]=deal(.5*(Ep(:,:,1)-(i/ep(3))*Ep(:,:,2)./khi),.5*((i*ep(3))*Ep(:,:,1).*khi+Ep(:,:,2)));
Ep(:,:,3)=-(i/ep(2))*repmat(u,ny,1).*Ep(:,:,1);
for ie=1:ne;
Ep(:,:,ie)=exp(-i*khi.*repmat(y(:),1,nu)).*Ep(:,:,ie);
Em(:,:,ie)=exp(i*khi.*repmat(y(:),1,nu)).*Em(:,:,ie);
end;
Ep(:,:,2:3)=-i*Ep(:,:,2:3);Em(:,:,2:3)=-i*Em(:,:,2:3);% de_clonnage


if cal_angles;
%   mise en forme des ondes planes  en 1D
km=[u(:),-khi(1,:).'];km=km./repmat(sqrt(sum((km.^2),2)),1,size(km,2));% vecteurs d'onde normalises
um=[km(:,2),-km(:,1)];% um= -ez vect.  km
HHm=repmat(um(:,1).',ny,1).*Em(:,:,2)+repmat(um(:,2).',ny,1).*Em(:,:,3);
kp=km;kp(:,2)=-kp(:,2);
up=um;up(:,1)=-up(:,1);
HHp=repmat(up(:,1).',ny,1).*Ep(:,:,2)+repmat(up(:,2).',ny,1).*Ep(:,:,3);
teta=atan2(real(kp(:,1)),real(kp(:,2)));
EEp=Ep(:,:,1).';EEm=Em(:,:,1).';% on met y en dernier
HHp=HHp.';HHm=HHm.';
Fp=pi*sqrt(ep(1)*ep(2))*retdiag(cos(teta).^2)*real(EEp.*conj(HHp));Fm=pi*sqrt(ep(1)*ep(2))*retdiag(cos(teta).^2)*real(EEm.*conj(HHm));
angles=struct('teta',teta,'kp',kp,'up',up,'EEp',EEp,'HHp',HHp,'Fp',full(Fp),'km',km,'um',um,'EEm',EEm,'HHm',HHm,'Fm',full(Fm));

end;
Ep=permute(Ep,[2,3,1]);Em=permute(Em,[2,3,1]);% on met y en dernier

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function varargout=ondes_planes_2D(n,u,v,e,x,y,z,wx,wy,wz,k0,parm);
if nargin<12;parm=[];end;
if nargin<11;k0=1;end;
if nargin<10;wz=[];end;
ep=ret2ep(n,k0);

if iscell(ep);ep=ret2ep(ep{:});end;

defaultopt=struct('apod',0,'delt',0,'uvmesh',0,'clone',0);
if nargin<11;parm=[];end;
if isempty(parm);parm=defaultopt;end;
apod=retoptimget(parm,'apod',defaultopt,'fast');
cal_angles=nargout>2;
delt=retoptimget(parm,'delt',defaultopt,'fast');
uvmesh=retoptimget(parm,'uvmesh',defaultopt,'fast');
clone=retoptimget(parm,'clone',defaultopt,'fast');
if clone==0;e(:,:,:,4:6)=i*e(:,:,:,4:6);end% clonage


if isempty(wx);% [ x y z ] -->[ y z x ]
e=permute(e,[2,3,1,4]);e=e(:,:,:,[2,3,1,5,6,4]);
[varargout{1:nargout}]=ondes_planes_2Dz(ep([2,3,1,5,6,4]),u,v,e,y,z,x,wy,wz,apod,cal_angles,delt,uvmesh);
if uvmesh==1;for ii=1:min(2,nargout);varargout{ii}=varargout{ii}(:,[3,1,2,6,4,5],:);end;% remise en ordre des composantes
else;for ii=1:min(2,nargout);varargout{ii}=varargout{ii}(:,:,[3,1,2,6,4,5],:);end;end;

return;end;
if isempty(wy);% [ x y z ] -->[ z x y]
e=permute(e,[3,1,2,4]);e=e(:,:,:,[3,1,2,6,4,5]);
[varargout{1:nargout}]=ondes_planes_2Dz(ep([3,1,2,6,4,5]),u,v,e,z,x,y,wz,wx,apod,cal_angles,delt,uvmesh);
if uvmesh==1;for ii=1:min(2,nargout);varargout{ii}=varargout{ii}(:,[2,3,1,5,6,4],:);end;% remise en ordre des composantes
else;for ii=1:min(2,nargout);varargout{ii}=varargout{ii}(:,:,[2,3,1,5,6,4],:);end;end;

return;end;
if isempty(wz);
[varargout{1:nargout}]=ondes_planes_2Dz(ep,u,v,e,x,y,z,wx,wy,apod,cal_angles,delt,uvmesh);
return;end;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [Ep,Em,angles]=ondes_planes_2Dz(ep,u,v,e,x,y,z,wx,wy,apod,cal_angles,delt,uvmesh);
  
nx=length(x);ny=length(y);nz=length(z);ne=size(e,4);
if imag(uvmesh)==0;% ****** pas cylindres popov
if apod~=0;
apodise_x=interp1(linspace(min(x),max(x),100),retchamp([apod,100]),x);
apodise_y=interp1(linspace(min(y),max(y),100),retchamp([apod,100]),y);
else;apodise_x=1;apodise_y=1;end;	
	
if uvmesh==1; % <----u v meshed
nuv=length(u);
Ep=zeros(nz,nuv,ne);Em=zeros(nz,nuv,ne);
e=retreshape(e,[nz,nx*ny,ne]);
wx=wx.*apodise_x/(2*pi);wy=wy.*apodise_y/(2*pi);
%wx=wx.*retchamp([apod,nx])/(2*pi);wy=wy.*retchamp([apod,ny])/(2*pi);
for iuv=1:nuv;%prv=ax(:,iuv)*ay(iuv,:);
prv=(exp(-i*u(iuv)*x(:)).*wx(:))*(exp(-i*v(iuv)*y(:)).*wy(:)).';
for ie=[1,2,4,5];
Ep(:,iuv,ie)=e(:,:,ie)*prv(:);
end;end;

clear e;
else;         % <----u v non  meshed
% ax=(exp(-i*u(:)*x(:).')*retdiag(wx.*retchamp([apod,nx]))/(2*pi)).';
% ay=exp(-i*v(:)*y(:).')*retdiag(wy.*retchamp([apod,ny]))/(2*pi);
ax=(exp(-i*u(:)*x(:).')*retdiag(wx.*apodise_x)/(2*pi)).';
ay=exp(-i*v(:)*y(:).')*retdiag(wy.*apodise_y)/(2*pi);
nu=length(u);nv=length(v);nuv=nu*nv;
Ep=zeros(nz,nu,nv,ne);Em=zeros(nz,nu,nv,ne);
for iu=1:nu;for iv=1:nv;
prv=ax(:,iu)*ay(iv,:);
prv=retreshape(repmat(full(prv(:).'),nz,1),[nz,nx,ny,1]);% full sinon bug si nx=1
for ie=[1,2,4,5];
Ep(:,iu,iv,ie)=sum(sum(e(:,:,:,ie).*prv,2),3);
end;
end;end;
clear e;
[u,v]=ndgrid(u,v);
Ep=retreshape(Ep,[nz,nuv,ne]);Em=retreshape(Em,[nz,nuv,ne]);
end;         % <----u v meshed ? 
else;                % ********** cylindres popov
uvmesh=real(uvmesh);
Ep=e;clear e;
if uvmesh==1; % <----u v meshed
nuv=length(u);
else;         % <----u v non  meshed
nu=length(u);nv=length(v);nuv=nu*nv;
[u,v]=ndgrid(u,v);
end;
end;

khi=retbidouille(retsqrt(ep(6)*ep(3)-u.^2-v.^2,-1));
khi=repmat(khi(:).',nz,1);
u=repmat(u(:).',nz,1);
v=repmat(v(:).',nz,1);

Ep(:,:,3)=(i/ep(6))*(u.*Ep(:,:,5)-v.*Ep(:,:,4));
Ep(:,:,6)=(i/ep(3))*(u.*Ep(:,:,2)-v.*Ep(:,:,1));

[Em(:,:,1),Em(:,:,2),Em(:,:,4),Em(:,:,5)]=deal(...
.5*( Ep(:,:,1)+i*(ep(2)*Ep(:,:,5)+i*u.*Ep(:,:,3))./khi),...
.5*( Ep(:,:,2)+i*(-ep(1)*Ep(:,:,4)+i*v.*Ep(:,:,3))./khi),...
.5*( Ep(:,:,4)+i*(ep(5)*Ep(:,:,2)+i*u.*Ep(:,:,6))./khi),...
.5*( Ep(:,:,5)+i*(-ep(4)*Ep(:,:,1)+i*v.*Ep(:,:,6))./khi));
Em(:,:,3)=(i/ep(6))*(u.*Em(:,:,5)-v.*Em(:,:,4));
Em(:,:,6)=(i/ep(3))*(u.*Em(:,:,2)-v.*Em(:,:,1));

[Ep(:,:,1),Ep(:,:,2),Ep(:,:,4),Ep(:,:,5)]=deal(...
.5*( Ep(:,:,1)-i*(ep(2)*Ep(:,:,5)+i*u.*Ep(:,:,3))./khi),...
.5*( Ep(:,:,2)-i*(-ep(1)*Ep(:,:,4)+i*v.*Ep(:,:,3))./khi),...
.5*( Ep(:,:,4)-i*(ep(5)*Ep(:,:,2)+i*u.*Ep(:,:,6))./khi),...
.5*( Ep(:,:,5)-i*(-ep(4)*Ep(:,:,1)+i*v.*Ep(:,:,6))./khi));
Ep(:,:,3)=(i/ep(6))*(u.*Ep(:,:,5)-v.*Ep(:,:,4));
Ep(:,:,6)=(i/ep(3))*(u.*Ep(:,:,2)-v.*Ep(:,:,1));

Ep(:,:,4:6)=-i*Ep(:,:,4:6);Em(:,:,4:6)=-i*Em(:,:,4:6);% de_clonnage
for ie=1:ne;
Ep(:,:,ie)=exp(-i*khi.*repmat(z(:),[1,nuv])).*Ep(:,:,ie);
Em(:,:,ie)=exp(i*khi.*repmat(z(:),[1,nuv])).*Em(:,:,ie);
end;
if cal_angles==1;
khi=khi(1,:);u=u(1,:);v=v(1,:);

km=[u(:),v(:),-khi(:)];km=km./repmat(sqrt(sum((km.^2),2)),1,size(km,2));%vecteurs d'onde normalises
vm=[-km(:,2),km(:,1),zeros(size(km,1),1)];
ii=find(abs(km(:,1).^2+km(:,2).^2)<100*eps);vm(ii,1)=-sin(delt);vm(ii,2)=cos(delt);%incidence normale: delt valeur par defaut de delt
vm=vm./repmat(sqrt(sum((vm.^2),2)),1,size(vm,2));
um=[km(:,3).*vm(:,2),-km(:,3).*vm(:,1),km(:,2).*vm(:,1)-km(:,1).*vm(:,2)];%um= vm vect.  km
kp=km;kp(:,3)=-kp(:,3);
vp=vm;up=um;up(:,1:2)=-up(:,1:2);

psip=zeros(nz,nuv);psim=zeros(nz,nuv);
EEp=zeros(nz,nuv,2);EEEp=zeros(nz,nuv,2);HHp=zeros(nz,nuv,2);HHHp=zeros(nz,nuv,2);
EEm=zeros(nz,nuv,2);EEEm=zeros(nz,nuv,2);HHm=zeros(nz,nuv,2);HHHm=zeros(nz,nuv,2);

%EEEEm=zeros(nz,nuv);HHHHm=zeros(nz,nuv);
vp=vp.';up=up.';kp=kp.';vm=vm.';um=um.';km=km.';

for iz=1:nz;

EEp(iz,:,1)=Ep(iz,:,1).*up(1,:)+Ep(iz,:,2).*up(2,:)+Ep(iz,:,3).*up(3,:);
HHp(iz,:,1)=Ep(iz,:,4).*up(1,:)+Ep(iz,:,5).*up(2,:)+Ep(iz,:,6).*up(3,:);
EEm(iz,:,1)=Em(iz,:,1).*um(1,:)+Em(iz,:,2).*um(2,:)+Em(iz,:,3).*um(3,:);
HHm(iz,:,1)=Em(iz,:,4).*um(1,:)+Em(iz,:,5).*um(2,:)+Em(iz,:,6).*um(3,:);

EEp(iz,:,2)=Ep(iz,:,1).*vp(1,:)+Ep(iz,:,2).*vp(2,:)+Ep(iz,:,3).*vp(3,:);
HHp(iz,:,2)=Ep(iz,:,4).*vp(1,:)+Ep(iz,:,5).*vp(2,:)+Ep(iz,:,6).*vp(3,:);
EEm(iz,:,2)=Em(iz,:,1).*vm(1,:)+Em(iz,:,2).*vm(2,:)+Em(iz,:,3).*vm(3,:);
HHm(iz,:,2)=Em(iz,:,4).*vm(1,:)+Em(iz,:,5).*vm(2,:)+Em(iz,:,6).*vm(3,:);

% analyse des polarisations
psip(iz,:)=.5*atan2(real(EEp(iz,:,1)).*real(EEp(iz,:,2))+imag(EEp(iz,:,1)).*imag(EEp(iz,:,2)),.5*(abs(EEp(iz,:,1)).^2-abs(EEp(iz,:,2)).^2));
EEEp(iz,:,1)=EEp(iz,:,1).*cos(psip(iz,:))+EEp(iz,:,2).*sin(psip(iz,:));EEEp(iz,:,2)=-EEp(iz,:,1).*sin(psip(iz,:))+EEp(iz,:,2).*cos(psip(iz,:));
ii=find(abs(EEEp(iz,:,2))>abs(EEEp(iz,:,1)));psip(iz,ii)=psip(iz,ii)+pi/2;[EEEp(iz,ii,1),EEEp(iz,ii,2)]=deal(EEEp(iz,ii,2),-EEEp(iz,ii,1));
ii=find(psip(iz,:)<=-pi/2+100*eps);psip(iz,ii)=psip(iz,ii)+pi;EEEp(iz,ii,:)=-EEEp(iz,ii,:);
HHHp(iz,:,1)=HHp(iz,:,1).*cos(psip(iz,:))+HHp(iz,:,2).*sin(psip(iz,:));HHHp(iz,:,2)=-HHp(iz,:,1).*sin(psip(iz,:))+HHp(iz,:,2).*cos(psip(iz,:));
ii=find(abs(psip(iz,:)-pi/2)<100*eps);psip(iz,ii)=pi/2;ii=find(abs(psip(iz,:))<100*eps);psip(iz,ii)=0;

psim(iz,:)=.5*atan2(real(EEm(iz,:,1)).*real(EEm(iz,:,2))+imag(EEm(iz,:,1)).*imag(EEm(iz,:,2)),.5*(abs(EEm(iz,:,1)).^2-abs(EEm(iz,:,2)).^2));
EEEm(iz,:,1)=EEm(iz,:,1).*cos(psim(iz,:))+EEm(iz,:,2).*sin(psim(iz,:));EEEm(iz,:,2)=-EEm(iz,:,1).*sin(psim(iz,:))+EEm(iz,:,2).*cos(psim(iz,:));
ii=find(abs(EEEm(iz,:,2))>abs(EEEm(iz,:,1)));psim(iz,ii)=psim(iz,ii)+pi/2;[EEEm(iz,ii,1),EEEm(iz,ii,2)]=deal(EEEm(iz,ii,2),-EEEm(iz,ii,1));
ii=find(psim(iz,:)<=-pi/2+100*eps);psim(iz,ii)=psim(iz,ii)+pi;EEEm(iz,ii,:)=-EEEm(iz,ii,:);
HHHm(iz,:,1)=HHm(iz,:,1).*cos(psim(iz,:))+HHm(iz,:,2).*sin(psim(iz,:));HHHm(iz,:,2)=-HHm(iz,:,1).*sin(psim(iz,:))+HHm(iz,:,2).*cos(psim(iz,:));
ii=find(abs(psim(iz,:)-pi/2)<100*eps);psim(iz,ii)=pi/2;ii=find(abs(psim(iz,:))<100*eps);psim(iz,ii)=0;

end; %iz

vp=vp.';up=up.';kp=kp.';vm=vm.';um=um.';km=km.';
teta=acos(kp(:,3));delta=atan2(-real(vm(:,1)),real(vm(:,2)));

if uvmesh~=1; % <----u v pas meshed
Ep=retreshape(Ep,nz,nu,nv,6);Em=retreshape(Em,nz,nu,nv,6);
EEp=retreshape(EEp,nz,nu,nv,2);EEEp=retreshape(EEEp,nz,nu,nv,2);HHp=retreshape(HHp,nz,nu,nv,2);HHHp=retreshape(HHHp,nz,nu,nv,2);
EEm=retreshape(EEm,nz,nu,nv,2);EEEm=retreshape(EEEm,nz,nu,nv,2);HHm=retreshape(HHm,nz,nu,nv,2);HHHm=retreshape(HHHm,nz,nu,nv,2);
psip=retreshape(psip,nz,nu,nv);kp=retreshape(kp,nu,nv,3);up=retreshape(up,nu,nv,3);vp=retreshape(vp,nu,nv,3);
psim=retreshape(psim,nz,nu,nv);km=retreshape(km,nu,nv,3);um=retreshape(um,nu,nv,3);vm=retreshape(vm,nu,nv,3);
% on met z en dernier
EEp=permute(EEp,[2,3,4,1]);EEm=permute(EEm,[2,3,4,1]);EEEp=permute(EEEp,[2,3,4,1]);EEEm=permute(EEEm,[2,3,4,1]);
HHp=permute(HHp,[2,3,4,1]);HHm=permute(HHm,[2,3,4,1]);HHHp=permute(HHHp,[2,3,4,1]);HHHm=permute(HHHm,[2,3,4,1]);
psip=permute(psip,[2,3,1]);psim=permute(psim,[2,3,1]);
Fp=2*pi^2*ep(1)*ep(4)*real(EEp(:,:,1,:).*conj(HHp(:,:,2,:))-EEp(:,:,2,:).*conj(HHp(:,:,1,:)));
Fm=2*pi^2*ep(1)*ep(4)*real(EEm(:,:,1,:).*conj(HHm(:,:,2,:))-EEm(:,:,2,:).*conj(HHm(:,:,1,:)));
Fp=retdiag(cos(teta).^2)*retreshape(Fp,[nuv,nz]);Fp=retreshape(Fp,[nu,nv,nz]);
Fm=retdiag(cos(teta).^2)*retreshape(Fm,[nuv,nz]);Fm=retreshape(Fm,[nu,nv,nz]);
teta=retreshape(teta,nu,nv);delta=retreshape(delta,nu,nv);


else;%   <----u v meshed
% on met z en dernier
EEp=permute(EEp,[2,3,1]);EEm=permute(EEm,[2,3,1]);EEEp=permute(EEEp,[2,3,1]);EEEm=permute(EEEm,[2,3,1]);
HHp=permute(HHp,[2,3,1]);HHm=permute(HHm,[2,3,1]);HHHp=permute(HHHp,[2,3,1]);HHHm=permute(HHHm,[2,3,1]);
psip=permute(psip,[2,1]);psim=permute(psim,[2,1]);

Fp=2*pi^2*ep(1)*ep(4)*real(EEp(:,1,:).*conj(HHp(:,2,:))-EEp(:,2,:).*conj(HHp(:,1,:)));Fp=retdiag(cos(teta).^2)*retreshape(Fp,[nuv,nz]);
Fm=2*pi^2*ep(1)*ep(4)*real(EEm(:,1,:).*conj(HHm(:,2,:))-EEm(:,2,:).*conj(HHm(:,1,:)));Fm=retdiag(cos(teta).^2)*retreshape(Fm,[nuv,nz]);

end;% <----u v meshed %

%angles={teta,delta,psip,kp,up,vp,EEp,HHp,EEEp,HHHp,psim,km,um,vm,EEm,HHm,EEEm,HHHm};
angles=struct('teta',teta,'delta',delta,'psip',psip,'kp',kp,'vp',vp,'up',up,'EEp',EEp,'HHp',HHp,'EEEp',EEEp,'HHHp',HHHp,'Fp',full(Fp),'psim',psim,'km',km,'vm',vm,'um',um,'EEm',EEm,'HHm',HHm,'EEEm',EEEm,'HHHm',HHHm,'Fm',full(Fm));
end;  % cal_angles
if uvmesh==1; % <----u v meshed
Ep=permute(Ep,[2,3,1]);Em=permute(Em,[2,3,1]);% on met z en dernier
else;
Ep=retreshape(Ep,nz,nu,nv,6);Em=retreshape(Em,nz,nu,nv,6);
Ep=permute(Ep,[2,3,4,1]);Em=permute(Em,[2,3,4,1]);% on met z en dernier
end;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [eep,eem]=ondes_planes_1D_ancien(u,ep,e,yy,ww);% ancienne version 1 D
%developpement en ondes planes dans un milieu homogene ep=[epsilon;mux;1/muy] 
% constante de propagation sur y:u 
%eep: x croissant eem x decroissant
%ez=e(:,1) hy=e(:,2) champs donnes aux points yy (poids d'integration ww) vecteurs ligne
eez=ww*((e(:,1)*ones(size(u))).*exp(-i.*yy.'*u));
hhy=ww*((e(:,3)*ones(size(u))).*exp(-i.*yy.'*u));
hhy=-ep(2)*hhy./sqrt(ep(1)*ep(2)-u.^2);
eep=(eez+hhy)./(4*pi);eem=(eez-hhy)./(4*pi);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function varargout=ondes_planes_popov(init,sh,sb,a,n,u,v,z,inc,k0,parm);
%  function [Ep,Em,angles]=retop(init,sh,sb,a,n,u,v,z,inc,k0,parm);
% modif Alex: ajout de inc


% if nargin < 11;
% 	% assert(nargin == 10);% Pas valable en Matlab 7.3
% 	parm = k0;
% 	k0 = inc;
% 	inc = 1;
% end
if nargin==9;[inc,k0,parm]=deal(1,inc,[]);end;% modif 3 2011
if nargin==10;if isstruct(k0)|isempty(k0);[inc,k0,parm]=deal(1,inc,k0);else;parm=[];end;end;


sh=retio(sh);sb=retio(sb);a=retio(a);

N=init{end}.nhanckel;sym=init{end}.sym;L=init{2};k=init{3};wk=init{4};cao=init{5};
ep_mu=ret2ep(n,k0);
if iscell(ep_mu);ep_mu=ret2ep(ep_mu{:});end;
defaultopt=struct('uvmesh',0,'delt',0);
%if nargin<10;parm=[];end;
if isempty(parm);parm=defaultopt;end;

uvmesh=retoptimget(parm,'uvmesh',defaultopt,'fast');
delt=retoptimget(parm,'delt',defaultopt,'fast');
if uvmesh==0;[uu,vv]=ndgrid(u,v);uu=uu(:);vv=vv(:);else;uu=u(:);vv=v(:);end;

kk=sqrt(uu.^2+vv.^2);phi_p=atan2(vv,uu)-pi/2;

Ep=zeros(length(z),numel(uu),6);
if sym==0;
expp=exp(i*(L+1)*phi_p);expm=exp(i*(L-1)*phi_p);
else;
cosp=cos((L+1)*phi_p);cosm=cos((L-1)*phi_p);sinp=sin((L+1)*phi_p);sinm=sin((L-1)*phi_p);
end;

if isfield(parm,'rmax');% calcul de champ sur un rayon puis transformation de Hankel 
Pml=init{10};
r=real(retinterp_popov(a{7}{1},Pml,1));
if ~isempty(Pml);xx=retelimine([parm.rmax,real(retinterp_popov(r,Pml,1)),0]);else;xx=retelimine([rmax,r,0]);end;% on ajoute les pml reelles aux discontinuites
xx=xx(xx<=parm.rmax*(1+10*eps));
[x,wx]=retgauss(0,max(xx),50,2,xx);
if z(1)==0;tab=[0,1,1];else;tab=[[0,1,1];[z(1),1,0]];end;for ii=1:length(z)-1;tab=[[0,1,1];[z(ii+1)-z(ii),1,0];tab];end;
% modif Alex: ajout de inc
iinit=init;iinit{end}.sym=0;e=retchamp(iinit,{a},sh,sb,inc,{x,0},tab);
eep=.5*(e(:,:,:,2)-i*e(:,:,:,1));eem=.5*(e(:,:,:,2)+i*e(:,:,:,1));
hhp=.5i*(e(:,:,:,5)-i*e(:,:,:,4));hhm=.5i*(e(:,:,:,5)+i*e(:,:,:,4));clear e;% clonage
JLp1=retbessel('j',L+1,kk(:)*x(:).')*retdiag(x.*wx);
JLm1=retbessel('j',L-1,kk(:)*x(:).')*retdiag(x.*wx);
end;

for kz=1:length(z);% <************** kz
if isfield(parm,'rmax');
ep=JLp1*retcolonne(eep(kz,:,:));em=JLm1*retcolonne(eem(kz,:,:));	
hp=JLp1*retcolonne(hhp(kz,:,:));hm=JLm1*retcolonne(hhm(kz,:,:));	
else;
% modif Alex: ajout de inc
uv=retsc(sh,retss(retc(a,z(kz)),sb),inc,2*N);
ep=retinterp(k,uv(1:N),kk,'pchip');em=retinterp(k,uv(N+1:2*N),kk,'pchip');
hp=retinterp(k,uv(2*N+1:3*N),kk,'pchip');hm=retinterp(k,uv(3*N+1:4*N),kk,'pchip');
end;
switch sym;% sym
case 0
Ep(kz,:,1)=(i/(2*pi))*(expp.*ep-expm.*em).';
Ep(kz,:,2)=(1/(2*pi))*(expp.*ep+expm.*em).';
Ep(kz,:,4)=(i/(2*pi))*(expp.*hp-expm.*hm).';
Ep(kz,:,5)=(1/(2*pi))*(expp.*hp+expm.*hm).';
case 1
Ep(kz,:,1)=(i/(2*pi))*(cosp.*ep-cosm.*em).';
Ep(kz,:,2)=(i/(2*pi))*(sinp.*ep+sinm.*em).';
Ep(kz,:,4)=(1/(2*pi))*(-sinp.*hp+sinm.*hm).';
Ep(kz,:,5)=(1/(2*pi))*(cosp.*hp+cosm.*hm).';
case -1
Ep(kz,:,1)=(1/(2*pi))*(-sinp.*ep+sinm.*em).';
Ep(kz,:,2)=(1/(2*pi))*(cosp.*ep+cosm.*em).';
Ep(kz,:,4)=(i/(2*pi))*(cosp.*hp-cosm.*hm).';
Ep(kz,:,5)=(i/(2*pi))*(sinp.*hp+sinm.*hm).';
end;% sym
end;% <************** kz
cal_angles=nargout>2;
[varargout{1:nargout}]=ondes_planes_2Dz(ep_mu,u,v,Ep,0,0,z,0,0,0,cal_angles,delt,uvmesh+i);






%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%             VERSIONS 2010                      %
%         integration locale                     %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%
%%%%%%  1D  %%%%%%%
%%%%%%%%%%%%%%%%%%%
function angles=caldop_1D(n,y,u,  e_x,x_x,y_x,w_x,   e_y,x_y,y_y,w_y, pol,k0,sens,parm);% parm ajouté 10 2011
if nargin<15;parm=[];end;
if nargin==14 & isstruct(sens);parm=sens;end;
if nargin<14 | (nargin==14 & (isstruct(sens)|isempty(sens)));[e_x,x_x,y_x,w_x,   e_y,x_y,y_y,w_y, pol,k0,sens]=deal(u,  e_x,x_x,y_x,w_x,   e_y,x_y,y_y,w_y, pol,k0);end;% si pas u (modes en 1D)
if imag(sens)~=0;%[ x y ] -->[ y x]
e_x=permute(e_x,[2,1,3]);e_x=e_x(:,:,[1,3,2]);e_x(:,:,1)=-e_x(:,:,1);
e_y=permute(e_y,[2,1,3]);e_y=e_y(:,:,[1,3,2]);e_y(:,:,1)=-e_y(:,:,1);
angles=caldop_1Dy(n,y,u,  e_y,y_y,x_y,w_y,  e_x,y_x,x_x,w_x,    pol,k0,real(sens),parm);
else;
angles=caldop_1Dy(n,y,u,  e_x,x_x,y_x,w_x,   e_y,x_y,y_y,w_y, pol,k0,sens,parm);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function angles=caldop_1Dy(n,y,u,  e_x,x_x,y_x,w_x,   e_y,x_y,y_y,w_y, pol,k0,sens,parm);
if length(n)>length(y)+1;
nb_couches=length(y)+1;
angles=caldop_modes_1D(n(nb_couches+1),n(1:nb_couches),y,u,  e_x,x_x,y_x,w_x,   e_y,x_y,y_y,w_y, pol,k0,sens,parm);
return;
end;
Flux_poynting=diff(retpoynting(e_x,[0,1],w_x,[]))+diff(retpoynting(e_y,[1,0],[],w_y));
% clonage des champs
e_x(:,:,2:3)=1i*e_x(:,:,2:3);
e_y(:,:,2:3)=1i*e_y(:,:,2:3);

init={pol,-u,k0};
sh=retb(init,n(1),1);sb=retb(init,n(end),-1);

Y=[y_x(:);y_y(:)];% tous les Y où on veut le champ
[Y_reduit,prv,kkY]=retelimine(Y,-1);
if isempty(y);y=mean(Y_reduit);n=[n(1),n(1)];end;%hh=max(Y_reduit)-max(y)+1.e-6;hb=min(y)-min(Y_reduit)+1.e-6;

hh=max(max(Y_reduit)-max(y),0)+1.e-6;hb=max(min(y)-min(Y_reduit),0)+1.e-6;% modif 6 2014



tab=[[hh;-diff(y(:));hb],n(:)];
if sens==1;inc=[0,1];n_op=n(1);y_op=y(1);else;inc=[1,0];n_op=n(end);y_op=y(end);end;

tab0=tab;tab0(:,2)=n_op;sh0=retb(init,n_op,1);sb0=retb(init,n_op,-1);
einc=retchamp(init,tab0,sh0,sb0,inc,{y_op-y(end)+hb});einc(:,:,2:3,:,:)=1i*einc(:,:,2:3,:,:);% clonage

e=retchamp(init,tab,sh,sb,inc,{Y_reduit-y(end)+hb});e(:,:,2:3,:,:)=1i*e(:,:,2:3,:,:);% clonage
for ii=1:size(e,5);e(:,:,:,:,ii)=e(:,:,:,:,ii)/einc(:,:,1,:,ii);end;% normalisation E=1

e=e(kkY,:,:,:,:);
ee_x=e(1:length(y_x),:,:,:,:);
ee_y=e(1+length(y_x):end,:,:,:,:);
EE=zeros(length(u),1);
mu=k0*n_op^pol;khi=sqrt(max(eps,(k0*n_op)^2-u.^2));
for ii=1:length(u);%<<<<<<<< ii
eee_x=zeros(size(e_x));
prv=exp(-1i*u(ii)*retcolonne(x_x,1));
for jj=1:3;
eee_x(:,:,jj)=ee_x(:,:,jj,1,ii)*prv;
end
eee_y=zeros(size(e_y));
prv=exp(-1i*u(ii)*retcolonne(x_y,1));
for jj=1:3;
eee_y(:,:,jj)=ee_y(:,:,jj,1,ii)*prv;
end
EE(ii)=-mu/(4i*pi*khi(ii))*(diff(ret_int_lorentz(e_x,eee_x,w_x,[]))+diff(ret_int_lorentz(e_y,eee_y,[],w_y)));
end;              %<<<<<<<< ii

HH=(k0*n_op/mu)*EE;% decloné
%   mise en forme des ondes planes  en 1D

if sens==1;
K=[u(:),khi(:)];K=K./repmat(sqrt(sum((K.^2),2)),1,size(K,2));% vecteurs d'onde normalises
U=[K(:,2),-K(:,1)];
else;
K=[u(:),-khi(:)];K=K./repmat(sqrt(sum((K.^2),2)),1,size(K,2));% vecteurs d'onde normalises
U=[K(:,2),-K(:,1)];
end
teta=atan2(u(:),khi(:));

F=pi*(k0*n_op)*(cos(teta).^2).*real(EE.*conj(HH));
e=full(exp(-.25i*pi)*sqrt(2*pi*k0*n_op)*retdiag(cos(teta))*[EE,HH.*U(:,1),HH.*U(:,2)]);% ajout 12 2014
angles=struct('teta',teta,'k',K,'u',U,'e',e,'EE',EE,'HH',HH,'F',full(F),'origine',y_op,'Flux_poynting',Flux_poynting);
if isfield(parm,'test') & parm.test==1;
if isfield(parm,'Scale');Scale=parm.Scale;else;Scale=1;end;
trace_boite_1D( e_x,x_x*Scale,y_x*Scale,w_x,   e_y,x_y*Scale,y_y*Scale,w_y ,'');		
end;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%
%%%%%%  2D  %%%%%%%
%%%%%%%%%%%%%%%%%%%


function angles=caldop_2D(n,z,u,v,  e_xy,x_xy,y_xy,z_xy,wx_xy,wy_xy,    e_zx,x_zx,y_zx,z_zx,wz_zx,wx_zx,    e_yz,x_yz,y_yz,z_yz,wy_yz,wz_yz, k0,sens,parm);
if nargin<25;parm=[];end;
if imag(sens)>0;% [ x y z ] -->[ y z x ]  strates perpendiculaires à Ox
e_xy=permute(e_xy,[2,3,1,4]);e_xy=e_xy(:,:,:,[2,3,1,5,6,4]);
e_zx=permute(e_zx,[2,3,1,4]);e_zx=e_zx(:,:,:,[2,3,1,5,6,4]);
e_yz=permute(e_yz,[2,3,1,4]);e_yz=e_yz(:,:,:,[2,3,1,5,6,4]);
angles=caldop_2Dz(n,z,u,v,  e_yz,y_yz,z_yz,x_yz,wy_yz,wz_yz,    e_xy,y_xy,z_xy,x_xy,wx_xy,wy_xy,    e_zx,y_zx,z_zx,x_zx,wz_zx,wx_zx, k0,real(sens),parm);

return
end;

if imag(sens)<0;% [ x y z ] -->[ z x y]  strates perpendiculaires à Oy
e_xy=permute(e_xy,[3,1,2,4]);e_xy=e_xy(:,:,:,[3,1,2,6,4,5]);
e_zx=permute(e_zx,[3,1,2,4]);e_zx=e_zx(:,:,:,[3,1,2,6,4,5]);
e_yz=permute(e_yz,[3,1,2,4]);e_yz=e_yz(:,:,:,[3,1,2,6,4,5]);
angles=caldop_2Dz(n,z,u,v,  e_zx,z_zx,x_zx,y_zx,wz_zx,wx_zx,    e_yz,z_yz,x_yz,y_yz,wy_yz,wz_yz,    e_xy,z_xy,x_xy,y_xy,wx_xy,wy_xy, k0,real(sens),parm);
return
end;

angles=caldop_2Dz(n,z,u,v,  e_xy,x_xy,y_xy,z_xy,wx_xy,wy_xy,    e_zx,x_zx,y_zx,z_zx,wz_zx,wx_zx,    e_yz,x_yz,y_yz,z_yz,wy_yz,wz_yz, k0,sens,parm);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function angles=caldop_2Dz(n,z,u,v,  e_xy,x_xy,y_xy,z_xy,wx_xy,wy_xy,    e_zx,x_zx,y_zx,z_zx,wz_zx,wx_zx,    e_yz,x_yz,y_yz,z_yz,wy_yz,wz_yz, k0,sens,parm);
if nargin<25;parm=[];end;
Flux_poynting=diff(retpoynting(e_xy,[0,0,1],wx_xy,wy_xy,[]))+diff(retpoynting(e_zx,[0,1,0],wx_zx,[],wz_zx))+diff(retpoynting(e_yz,[1,0,0],[],wy_yz,wz_yz));

if size(n,4)>1;% mode de Choon
angles=caldop_modes_2D_choon(n,z,u,v,  e_xy,x_xy,y_xy,z_xy,wx_xy,wy_xy,    e_zx,x_zx,y_zx,z_zx,wz_zx,wx_zx,    e_yz,x_yz,y_yz,z_yz,wy_yz,wz_yz, k0,sens,parm);
return;
end;    

if length(n)>length(z)+1;
nb_couches=length(z)+1;
angles=caldop_modes_2D(n(nb_couches+1),n(nb_couches+2),n(1:nb_couches),z,u,v,  e_xy,x_xy,y_xy,z_xy,wx_xy,wy_xy,    e_zx,x_zx,y_zx,z_zx,wz_zx,wx_zx,    e_yz,x_yz,y_yz,z_yz,wy_yz,wz_yz, k0,sens,parm);
return;
end;
defaultopt=struct('uvmesh',0,'test',0,'option_parfor',0);
if isempty(parm);parm=defaultopt;end;
uvmesh=retoptimget(parm,'uvmesh',defaultopt,'fast');
test=retoptimget(parm,'test',defaultopt,'fast');
option_parfor=retoptimget(parm,'option_parfor',defaultopt,'fast');
if ~uvmesh;
nu=length(u);nv=length(v);nuv=nu*nv;
[u,v]=ndgrid(u,v);
end;

u=u(:);v=v(:);

if sens==1;inc=[0,1];n_op=n(1);z_op=z(1);else;inc=[1,0];n_op=n(end);z_op=z(end);end;
khi=sqrt(max(eps,(k0*n_op)^2-u.^2-v.^2));
delt=0;
	if sens==1;
	K=[u,v,khi];K=K./repmat(sqrt(sum((K.^2),2)),1,size(K,2));%vecteurs d'onde normalises
	V=[-K(:,2),K(:,1),zeros(size(K,1),1)];
	ii=find(abs(K(:,1).^2+K(:,2).^2)<100*eps);V(ii,1)=-sin(delt);V(ii,2)=cos(delt);%incidence normale: delt valeur par defaut de delt
	V=V./repmat(sqrt(sum((V.^2),2)),1,size(V,2));
	U=[K(:,3).*V(:,2),-K(:,3).*V(:,1),K(:,2).*V(:,1)-K(:,1).*V(:,2)];%U= V vect.  K
	else;
	K=[u,v,-khi];K=K./repmat(sqrt(sum((K.^2),2)),1,size(K,2));%vecteurs d'onde normalises
	V=[-K(:,2),K(:,1),zeros(size(K,1),1)];
	ii=find(abs(K(:,1).^2+K(:,2).^2)<100*eps);V(ii,1)=-sin(delt);V(ii,2)=cos(delt);%incidence normale: delt valeur par defaut de delt
	V=V./repmat(sqrt(sum((V.^2),2)),1,size(V,2));
	U=[K(:,3).*V(:,2),-K(:,3).*V(:,1),K(:,2).*V(:,1)-K(:,1).*V(:,2)];%U= V vect.  K
	end;
delta=atan2(-V(:,1),V(:,2));cos_delta=cos(delta);sin_delta=sin(delta);
[uv_reduit,prv,kkuv]=retelimine(sqrt(u.^2+v.^2),-1);

Z=[z_xy(:);z_zx(:);z_yz(:)];% tous les Z où on veut le champ
[Z_reduit,prv,kkZ]=retelimine(Z,-1);
if isempty(z);z=mean(Z_reduit);n=[n(1),n(1)];end;
%hh=max(Z_reduit)-z(1);hb=z(end)-min(Z_reduit);% modif 2011
%hh=max(Z_reduit)-max(z)+1.e-6;hb=min(z)-min(Z_reduit)+1.e-6;
hh=max(max(Z_reduit)-max(z),0)+1.e-6;hb=max(min(z)-min(Z_reduit),0)+1.e-6;% modif 6 2014


tab=[[hh;-diff(z(:));hb],n(:)];
tab0=tab;tab0(:,2)=n_op;

e=cell(1,2);
for kpol=1:2;
init={2*kpol-2,-uv_reduit,k0};
sh=retb(init,n(1),1);sb=retb(init,n(end),-1);
sh_inc=retb(init,n_op,1);sb_inc=retb(init,n_op,-1);
e_inc=retchamp(init,tab0,sh_inc,sb_inc,inc,{z_op-z(end)+hb});e_inc(:,:,2:3,:,:)=1i*e_inc(:,:,2:3,:,:);% clonage
e{kpol}=retchamp(init,tab,sh,sb,inc,{Z_reduit-z(end)+hb});e{kpol}(:,:,2:3,:,:)=1i*e{kpol}(:,:,2:3,:,:);% clonage
for ii=1:size(e{kpol},5);e{kpol}(:,:,:,:,ii)=e{kpol}(:,:,:,:,ii)/e_inc(:,:,1,:,ii);end;% normalisation E=1
end;  % kpol


% clonage des champs
e_xy(:,:,:,4:6)=1i*e_xy(:,:,:,4:6);
e_zx(:,:,:,4:6)=1i*e_zx(:,:,:,4:6);
e_yz(:,:,:,4:6)=1i*e_yz(:,:,:,4:6);

[X_xy,Y_xy]=ndgrid(x_xy,y_xy);
[X_zx,Y_zx]=ndgrid(x_zx,y_zx);
[X_yz,Y_yz]=ndgrid(x_yz,y_yz);
EE=zeros(length(u),2);
HH=zeros(length(u),2);

for kpol=1:2;
if option_parfor
parfor ii=1:length(u);
eee_xy=cale(u(ii),v(ii),X_xy,Y_xy,e{kpol}(kkZ(1:length(z_xy)),:,:,:,kkuv(ii)),size(e_xy),cos_delta(ii),sin_delta(ii),kpol);
eee_zx=cale(u(ii),v(ii),X_zx,Y_zx,e{kpol}(kkZ(1+length(z_xy):length(z_xy)+length(z_zx)),:,:,:,kkuv(ii)),size(e_zx),cos_delta(ii),sin_delta(ii),kpol);
eee_yz=cale(u(ii),v(ii),X_yz,Y_yz,e{kpol}(kkZ(1+length(z_xy)+length(z_zx):length(z_xy)+length(z_zx)+length(z_yz)),:,:,:,kkuv(ii)),size(e_yz),cos_delta(ii),sin_delta(ii),kpol);
if kpol==1;
EE(ii,2)=k0/(8i*pi^2*khi(ii))*(diff(ret_int_lorentz(e_xy,eee_xy,wx_xy,wy_xy,[]))+diff(ret_int_lorentz(e_yz,eee_yz,[],wy_yz,wz_yz))+diff(ret_int_lorentz(e_zx,eee_zx,wx_zx,[],wz_zx)));
else;
HH(ii,2)=k0*n_op^2/(8i*pi^2*khi(ii))*(diff(ret_int_lorentz(e_xy,eee_xy,wx_xy,wy_xy,[]))+diff(ret_int_lorentz(e_yz,eee_yz,[],wy_yz,wz_yz))+diff(ret_int_lorentz(e_zx,eee_zx,wx_zx,[],wz_zx)));
end;	
end;%ii
else;
    for ii=1:length(u);
eee_xy=cale(u(ii),v(ii),X_xy,Y_xy,e{kpol}(kkZ(1:length(z_xy)),:,:,:,kkuv(ii)),size(e_xy),cos_delta(ii),sin_delta(ii),kpol);
eee_zx=cale(u(ii),v(ii),X_zx,Y_zx,e{kpol}(kkZ(1+length(z_xy):length(z_xy)+length(z_zx)),:,:,:,kkuv(ii)),size(e_zx),cos_delta(ii),sin_delta(ii),kpol);
eee_yz=cale(u(ii),v(ii),X_yz,Y_yz,e{kpol}(kkZ(1+length(z_xy)+length(z_zx):length(z_xy)+length(z_zx)+length(z_yz)),:,:,:,kkuv(ii)),size(e_yz),cos_delta(ii),sin_delta(ii),kpol);
if kpol==1;
EE(ii,2)=k0/(8i*pi^2*khi(ii))*(diff(ret_int_lorentz(e_xy,eee_xy,wx_xy,wy_xy,[]))+diff(ret_int_lorentz(e_yz,eee_yz,[],wy_yz,wz_yz))+diff(ret_int_lorentz(e_zx,eee_zx,wx_zx,[],wz_zx)));
else;
HH(ii,2)=k0*n_op^2/(8i*pi^2*khi(ii))*(diff(ret_int_lorentz(e_xy,eee_xy,wx_xy,wy_xy,[]))+diff(ret_int_lorentz(e_yz,eee_yz,[],wy_yz,wz_yz))+diff(ret_int_lorentz(e_zx,eee_zx,wx_zx,[],wz_zx)));
end;	
end;%ii

    end;
end;%kpol

mu=k0;ep=k0*n_op^2;
HH(:,1)=-1i*(k0*n_op/mu)*EE(:,2);% -1 à cause de l'orientation
EE(:,1)=-1i*(k0*n_op/ep)*HH(:,2);

HH=-1i*HH;% declonage
% analyse des polarisations
psi=.5*atan2(real(EE(:,1)).*real(EE(:,2))+imag(EE(:,1)).*imag(EE(:,2)),.5*(abs(EE(:,1)).^2-abs(EE(:,2)).^2));
EEE(:,1)=EE(:,1).*cos(psi)+EE(:,2).*sin(psi);EEE(:,2)=-EE(:,1).*sin(psi)+EE(:,2).*cos(psi);
ii=find(abs(EEE(:,2))>abs(EEE(:,1)));psi(ii)=psi(ii)+pi/2;[EEE(ii,1),EEE(ii,2)]=deal(EEE(ii,2),-EEE(ii,1));
ii=find(psi<=-pi/2+100*eps);psi(ii)=psi(ii)+pi;EEE(ii,:)=-EEE(ii,:);
HHH(:,1)=HH(:,1).*cos(psi)+HH(:,2).*sin(psi);HHH(:,2)=-HH(:,1).*sin(psi(:))+HH(:,2).*cos(psi);
ii=find(abs(psi-pi/2)<100*eps);psi(ii)=pi/2;ii=find(abs(psi)<100*eps);psi(ii)=0;
teta=acos(K(:,3));

F=2*pi^2*(k0*n_op)^2*retdiag(cos(teta).^2)*real(EE(:,1).*conj(HH(:,2))-EE(:,2).*conj(HH(:,1)));
if ~uvmesh;
EE=retreshape(EE,nu,nv,[]);
HH=retreshape(HH,nu,nv,[]);
EEE=retreshape(EEE,nu,nv,[]);
HHH=retreshape(HHH,nu,nv,[]);
F=retreshape(F,nu,nv,[]);
end;
e=sens*2i*pi*n_op*k0*retdiag(cos(teta))* [ (retdiag(EE(:,1))*U+retdiag(EE(:,2))*V),(retdiag(HH(:,1))*U+retdiag(HH(:,2))*V)];
if sens>0;theta=teta;else;theta=pi-teta;end;
E=-2i*pi*n_op*k0*retdiag(cos(theta))* [ (retdiag(EE(:,1))*U+retdiag(EE(:,2))*V),(retdiag(HH(:,1))*U+retdiag(HH(:,2))*V)];

angles=struct('teta',teta,'theta',theta,'delta',delta,'psi',psi,'k',K,'v',V,'u',U,'e',e,'E',full(E),'EE',EE,'HH',HH,'EEE',EEE,'HHH',HHH,'F',full(F),'origine',z_op,'Flux_poynting',Flux_poynting);
%if isfield(parm,'test') & parm.test==1;
if test;
if isfield(parm,'Scale');Scale=parm.Scale;else;Scale=1;end;
trace_boite_2D(  e_xy,x_xy*Scale,y_xy*Scale,z_xy*Scale,    e_zx,x_zx*Scale,y_zx*Scale,z_zx*Scale,    e_yz,x_yz*Scale,y_yz*Scale,z_yz*Scale ,'');		
end;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function eee=cale1(u,v,X,Y,e,sz,cos_delta,sin_delta,kpol);
eee=zeros(sz);
prv=retcolonne(exp(-1i*(u*retcolonne(X,1)+v*retcolonne(Y,1))),1);% champ se propageant en -u -v
if kpol==1;
eee(:,:,:,2)=-retreshape(e(:,1,1)*prv,sz(1:3));
eee(:,:,:,4)=retreshape(e(:,1,2)*prv,sz(1:3));
eee(:,:,:,6)=retreshape(e(:,1,3)*prv,sz(1:3));
eee(:,:,:,1)=-eee(:,:,:,2)*sin_delta;eee(:,:,:,2)=eee(:,:,:,2)*cos_delta;
eee(:,:,:,5)=eee(:,:,:,4)*sin_delta;eee(:,:,:,4)=eee(:,:,:,4)*cos_delta;
else;
eee(:,:,:,5)=-retreshape(e(:,1,1)*prv,sz(1:3));
eee(:,:,:,1)=retreshape(e(:,1,2)*prv,sz(1:3));
eee(:,:,:,3)=retreshape(e(:,1,3)*prv,sz(1:3));
eee(:,:,:,2)=eee(:,:,:,1)*sin_delta;eee(:,:,:,1)=eee(:,:,:,1)*cos_delta;
eee(:,:,:,4)=-eee(:,:,:,5)*sin_delta;eee(:,:,:,5)=eee(:,:,:,5)*cos_delta;
end;
% plus rapide de le faire plus haut
% [eee(:,:,:,1),eee(:,:,:,2)]=deal(eee(:,:,:,1)*cos_delta-eee(:,:,:,2)*sin_delta,eee(:,:,:,1)*sin_delta+eee(:,:,:,2)*cos_delta);
% [eee(:,:,:,4),eee(:,:,:,5)]=deal(eee(:,:,:,4)*cos_delta-eee(:,:,:,5)*sin_delta,eee(:,:,:,4)*sin_delta+eee(:,:,:,5)*cos_delta);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function eee=cale(u,v,X,Y,e,sz,cos_delta,sin_delta,kpol);
e=retreshape(e,[],3);
prv=retcolonne(exp(-1i*(u*retcolonne(X,1)+v*retcolonne(Y,1))),1);% champ se propageant en -u -v
eee=zeros(size(e,1),size(prv,2),6);

if kpol==1;
eee(:,:,2)=-e(:,1)*prv;
eee(:,:,4)=e(:,2)*prv;
eee(:,:,6)=e(:,3)*prv;
eee(:,:,1)=-eee(:,:,2)*sin_delta;eee(:,:,2)=eee(:,:,2)*cos_delta;
eee(:,:,5)=eee(:,:,4)*sin_delta;eee(:,:,4)=eee(:,:,4)*cos_delta;
else;
eee(:,:,5)=-e(:,1)*prv;
eee(:,:,1)=e(:,2)*prv;
eee(:,:,3)=e(:,3)*prv;
eee(:,:,2)=eee(:,:,1)*sin_delta;eee(:,:,1)=eee(:,:,1)*cos_delta;
eee(:,:,4)=-eee(:,:,5)*sin_delta;eee(:,:,5)=eee(:,:,5)*cos_delta;
end;
eee=retreshape(eee,sz);

% plus rapide de le faire plus haut
% [eee(:,:,:,1),eee(:,:,:,2)]=deal(eee(:,:,:,1)*cos_delta-eee(:,:,:,2)*sin_delta,eee(:,:,:,1)*sin_delta+eee(:,:,:,2)*cos_delta);
% [eee(:,:,:,4),eee(:,:,:,5)]=deal(eee(:,:,:,4)*cos_delta-eee(:,:,:,5)*sin_delta,eee(:,:,:,4)*sin_delta+eee(:,:,:,5)*cos_delta);


%%%%%%%%%%%%%%%
%   MODES     %
%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%		
function trace_boite_2D(  e_xy,x_xy,y_xy,z_xy,    e_zx,x_zx,y_zx,z_zx,    e_yz,x_yz,y_yz,z_yz, texte);		
try;figure;hold on;
[X,Z]=ndgrid(x_zx,z_zx);ee=permute(sum(abs(e_zx(:,:,:,1:3).^2),4),[2,3,1]);surface(X,ones(size(Z))*y_zx(1),Z,squeeze(sqrt(ee(:,1,:))));surface(X,ones(size(Z))*y_zx(2),Z,squeeze(sqrt(ee(:,2,:))));
[Y,Z]=ndgrid(y_yz,z_yz);ee=permute(sum(abs(e_yz(:,:,:,1:3).^2),4),[2,3,1]);surface(ones(size(Z))*x_yz(1),Y,Z,squeeze(sqrt(ee(1,:,:))));surface(ones(size(Z))*x_yz(2),Y,Z,squeeze(sqrt(ee(2,:,:))));
[X,Y]=ndgrid(x_xy,y_xy);ee=permute(sum(abs(e_xy(:,:,:,1:3).^2),4),[2,3,1]);surface(X,Y,ones(size(X))*z_xy(1),squeeze(sqrt(ee(:,:,1))));surface(X,Y,ones(size(X))*z_xy(2),squeeze(sqrt(ee(:,:,2))));
plot3(x_yz([1,1,1,1,1]),y_zx([1,1,2,2,1]),z_xy([1,2,2,1,1]),'-g','linewidth',2);
plot3(x_yz([2,2,2,2,2]),y_zx([1,1,2,2,1]),z_xy([1,2,2,1,1]),'-g','linewidth',2);
plot3(x_yz([1,1,2,2,1]),y_zx([1,1,1,1,1]),z_xy([1,2,2,1,1]),'-g','linewidth',2);
plot3(x_yz([1,1,2,2,1]),y_zx([2,2,2,2,2]),z_xy([1,2,2,1,1]),'-g','linewidth',2);
shading flat;axis tight;axis equal;colormap hot;xlabel('x');ylabel('y');zlabel('z');view([-130,25]);
if nargin>12;title(texte);end;retfont(gcf,0);
catch;disp('Probleme dans le tracé du champ sur la boite');
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%		
function trace_boite_1D(   e_x,x_x,y_x,w_x,   e_y,x_y,y_y,w_y, texte);		
try;figure;hold on;
[X,Y]=ndgrid(x_x,y_x);
ee=abs(e_x(:,:,1).^2);plot3(X(:,1),Y(:,1),ee(1,:),'-r',X(:,2),Y(:,2),ee(2,:),'-r',X(:,1),Y(:,1),0*ee(1,:),'-k',X(:,2),Y(:,2),0*ee(2,:),'-k','linewidth',2);
ee=abs(e_x(:,:,2).^2);plot3(X(:,1),Y(:,1),ee(1,:),'--g',X(:,2),Y(:,2),ee(2,:),'--g','linewidth',2);
ee=abs(e_x(:,:,3).^2);plot3(X(:,1),Y(:,1),ee(1,:),':b',X(:,2),Y(:,2),ee(2,:),':b','linewidth',2);
[Y,X]=ndgrid(y_y,x_y);
ee=abs(e_y(:,:,1).^2);plot3(X(:,1),Y(:,1),ee(:,1),'-r',X(:,2),Y(:,2),ee(:,2),'-r',X(:,1),Y(:,1),0*ee(:,1),'-k',X(:,2),Y(:,2),0*ee(:,2),'-k','linewidth',2);
ee=abs(e_y(:,:,2).^2);plot3(X(:,1),Y(:,1),ee(:,1),'--g',X(:,2),Y(:,2),ee(:,2),'--g','linewidth',2);
ee=abs(e_y(:,:,3).^2);plot3(X(:,1),Y(:,1),ee(:,1),':b',X(:,2),Y(:,2),ee(:,2),':b','linewidth',2);
axis tight;xlabel('x');ylabel('y');view([-130,25]);
prvx=xlim;prvy=ylim;lx=diff(prvx);cx=mean(prvx);ly=diff(prvy);cy=mean(prvy);if lx>ly;ylim([cy-lx/2,cy+lx/2]);else;xlim([cx-ly/2,cx+ly/2]);end;
if nargin>8;title(texte);end;retfont(gcf,0);
catch;disp('Probleme dans le tracé du champ sur la boite');
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%		
function trace_boite_popov(   e_x,x_x,z_x,wx_x,   e_z,x_z,z_z,wz_z, texte);		
try;figure;hold on;
[X,Z]=ndgrid(x_x,z_x);
%ee=squeeze(sum(abs(e_x(:,:,:,1:3).^2),4));plot3(X(:,1),Z(:,1),ee(1,:),'-r',X(:,2),Z(:,2),ee(2,:),'-r',X(:,1),Z(:,1),0*ee(1,:),'-k',X(:,2),Z(:,2),0*ee(2,:),'-k','linewidth',2);
ee=squeeze(sum(abs(e_x(:,:,:,1).^2),4));plot3(X(:,1),Z(:,1),ee(1,:),'-r',X(:,2),Z(:,2),ee(2,:),'-r',X(:,1),Z(:,1),0*ee(1,:),'-k',X(:,2),Z(:,2),0*ee(2,:),'-k','linewidth',2);
ee=squeeze(sum(abs(e_x(:,:,:,2).^2),4));plot3(X(:,1),Z(:,1),ee(1,:),'--g',X(:,2),Z(:,2),ee(2,:),'--g','linewidth',2);
ee=squeeze(sum(abs(e_x(:,:,:,3).^2),4));plot3(X(:,1),Z(:,1),ee(1,:),':b',X(:,2),Z(:,2),ee(2,:),':b','linewidth',2);
[Z,X]=ndgrid(z_z,x_z);
%ee=squeeze(sum(abs(e_z(:,:,:,1:3).^2),4));plot3(X(:,1),Z(:,1),ee,'-r',X(:,1),Z(:,1),0*ee(:,1),'-k','linewidth',2);
ee=squeeze(sum(abs(e_z(:,:,:,1).^2),4));plot3(X(:,1),Z(:,1),ee,'-r',X(:,1),Z(:,1),0*ee(:,1),'-k','linewidth',2);
ee=squeeze(sum(abs(e_z(:,:,:,2).^2),4));plot3(X(:,1),Z(:,1),ee,'--g','linewidth',2);
ee=squeeze(sum(abs(e_z(:,:,:,3).^2),4));plot3(X(:,1),Z(:,1),ee,':b','linewidth',2);
axis tight;xlabel('x');ylabel('z');view([50,45]);
prvx=xlim;prvy=ylim;lx=diff(prvx);cx=mean(prvx);ly=diff(prvy);cy=mean(prvy);if lx>ly;ylim([cy-lx/2,cy+lx/2]);else;xlim([cx-ly/2,cx+ly/2]);end;
if nargin>8;title(texte);end;retfont(gcf,0);
catch;disp('Probleme dans le tracé du champ sur la boite');
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function angles=caldop_modes_1D(neff,n,y,u,  e_x,x_x,y_x,w_x,   e_y,x_y,y_y,w_y, pol,k0,sens,parm);
Flux_poynting=diff(retpoynting(e_x,[0,1],w_x,[]))+diff(retpoynting(e_y,[1,0],[],w_y));

% clonage des champs
e_x(:,:,2:3)=1i*e_x(:,:,2:3);
e_y(:,:,2:3)=1i*e_y(:,:,2:3);
[gama,ex,ey]=calmode(n,y,neff,pol,k0,y_x,y_y);
ee_x=expand_1D(ex,x_x,gama,-1);
ee_y=expand_1D(ey,x_y,gama,-1);
amp_p=.25i*(diff(ret_int_lorentz(e_x,ee_x,w_x,[]))+diff(ret_int_lorentz(e_y,ee_y,[],w_y)));
ee_x=expand_1D(ex,x_x,gama,1);
ee_y=expand_1D(ey,x_y,gama,1);
amp_m=.25i*(diff(ret_int_lorentz(e_x,ee_x,w_x,[]))+diff(ret_int_lorentz(e_y,ee_y,[],w_y)));

ex0=expand_1D(ex,0,gama,-1);ey0=expand_1D(ey,0,gama,-1);% mode à l'origine
ex0(:,:,2:3)=-i*ex0(:,:,2:3);ey0(:,:,2:3)=-i*ey0(:,:,2:3);% declonage
angles=struct('amp_p',amp_p,'Fp',abs(amp_p)^2,'amp_m',amp_m,'Fm',abs(amp_m)^2,'ex',ex0,'ey',ey0,'neff',gama/k0,'Flux_poynting',Flux_poynting);
if isfield(parm,'test') & parm.test==1;
if isfield(parm,'Scale');Scale=parm.Scale;else;Scale=1;end;
trace_boite_1D( e_x,x_x*Scale,y_x*Scale,w_x,   e_y,x_y*Scale,y_y*Scale,w_y ,'');		
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function angles=caldop_modes_2D_choon(e_mode,x_mode,z_mode,gama,  e_xy,x_xy,y_xy,z_xy,wx_xy,wy_xy,    e_zx,x_zx,y_zx,z_zx,wz_zx,wx_zx,    e_yz,x_yz,y_yz,z_yz,wy_yz,wz_yz, k0,sens,parm);
Flux_poynting=diff(retpoynting(e_xy,[0,0,1],wx_xy,wy_xy,[]))+diff(retpoynting(e_zx,[0,1,0],wx_zx,[],wz_zx))+diff(retpoynting(e_yz,[1,0,0],[],wy_yz,wz_yz));
gama=gama*k0;

% clonage des champs
e_xy(:,:,:,4:6)=1i*e_xy(:,:,:,4:6);
e_zx(:,:,:,4:6)=1i*e_zx(:,:,:,4:6);
e_yz(:,:,:,4:6)=1i*e_yz(:,:,:,4:6);
e_mode(:,:,:,4:6)=1i*e_mode(:,:,:,4:6);
[x_mode,iix]=retelimine(x_mode);[z_mode,iiz]=retelimine(z_mode);e_mode=e_mode(iiz,iix,1,:);
[ee_xy_p,ee_xy_m]=expand_mode(e_mode,x_mode,z_mode,x_xy,y_xy,z_xy,gama,sens);
[ee_zx_p,ee_zx_m]=expand_mode(e_mode,x_mode,z_mode,x_zx,y_zx,z_zx,gama,sens);
[ee_yz_p,ee_yz_m]=expand_mode(e_mode,x_mode,z_mode,x_yz,y_yz,z_yz,gama,sens);
amp_p=.25i*(diff(ret_int_lorentz(e_xy,ee_xy_p,wx_xy,wy_xy,[]))+diff(ret_int_lorentz(e_zx,ee_zx_p,wx_zx,[],wz_zx))+diff(ret_int_lorentz(e_yz,ee_yz_p,[],wy_yz,wz_yz)));
amp_m=.25i*(diff(ret_int_lorentz(e_xy,ee_xy_m,wx_xy,wy_xy,[]))+diff(ret_int_lorentz(e_zx,ee_zx_m,wx_zx,[],wz_zx))+diff(ret_int_lorentz(e_yz,ee_yz_m,[],wy_yz,wz_yz)));
angles=struct('amp_p',amp_p,'Fp',abs(amp_p)^2,'amp_m',amp_m,'Fm',abs(amp_m)^2,'Flux_poynting',Flux_poynting,'neff',gama/k0);
if isfield(parm,'test') & parm.test==1;
trace_boite_2D(  e_xy,x_xy,y_xy,z_xy,    e_zx,x_zx,y_zx,z_zx,    e_yz,x_yz,y_yz,z_yz ,'Champ total');		
trace_boite_2D(  ee_xy_p,x_xy,y_xy,z_xy,    ee_zx_p,x_zx,y_zx,z_zx,    ee_yz_p,x_yz,y_yz,z_yz,'Mode p interpolé');		
% trace_boite_2D(  ee_xy_m,x_xy,y_xy,z_xy,    ee_zx_m,x_zx,y_zx,z_zx,    ee_yz_m,x_yz,y_yz,z_yz,'Mode m interpolé');		
lz_xy=ret_int_lorentz(ee_xy_p,ee_xy_m,wx_xy,wy_xy,[]),lz_zx=ret_int_lorentz(ee_zx_p,ee_zx_m,wx_zx,[],wz_zx),lz_yz=ret_int_lorentz(ee_yz_p,ee_yz_m,[],wy_yz,wz_yz)
end;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [ee_p,ee_m]=expand_mode(e_mode,x_mode,z_mode,x,y,z,gama,sens);
if sens>0;% le mode se propage en y
y=y(:).';
[ee_p,ee_m]=deal(zeros(length(z)*length(x),length(y),6));
exp_p=exp(1i*gama*y);exp_m=exp(-1i*gama*y);
[Z,X]=ndgrid(z,x);
for ii=1:6;
eee=interp2(z_mode,x_mode,e_mode(:,:,1,ii).',Z,X);
ee_p(:,:,ii)=eee(:)*exp_p;
ee_m(:,:,ii)=eee(:)*exp_m;
end;
ee_p=retreshape(ee_p,length(z),length(x),length(y),6);
ee_m=retreshape(ee_m,length(z),length(x),length(y),6);
ee_m(:,:,:,[2,4,6])=-ee_m(:,:,:,[2,4,6]);

else;% le mode se propage en x
e_mode(:,:,:,[1,4])=-e_mode(:,:,:,[1,4]);y_mode=-x_mode;
e_mode=e_mode(:,:,:,[2,1,3,5,4,6]);
[ee_p,ee_m]=deal(zeros(length(z)*length(y),length(x),6));
exp_p=exp(1i*gama*x);exp_m=exp(-1i*gama*x);
[Z,Y]=ndgrid(z,y);
for ii=1:6;
eee=interp2(z_mode,y_mode,e_mode(:,:,1,ii).',Z,Y);
ee_p(:,:,ii)=eee(:)*exp_p;
ee_m(:,:,ii)=eee(:)*exp_m;
end;
ee_p=permute(retreshape(ee_p,length(z),length(y),length(x),6),[1,3,2,4]);
ee_m=permute(retreshape(ee_m,length(z),length(y),length(x),6),[1,3,2,4]);
ee_m(:,:,:,[1,5,6])=-ee_m(:,:,:,[1,5,6]);
end;% sens

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function angles=caldop_modes_2D(neff,pol,n,z,L,teta,  e_xy,x_xy,y_xy,z_xy,wx_xy,wy_xy,    e_zx,x_zx,y_zx,z_zx,wz_zx,wx_zx,    e_yz,x_yz,y_yz,z_yz,wy_yz,wz_yz, k0,sens,parm);
Flux_poynting=diff(retpoynting(e_xy,[0,0,1],wx_xy,wy_xy,[]))+diff(retpoynting(e_zx,[0,1,0],wx_zx,[],wz_zx))+diff(retpoynting(e_yz,[1,0,0],[],wy_yz,wz_yz));
% clonage des champs
e_xy(:,:,:,4:6)=1i*e_xy(:,:,:,4:6);
e_zx(:,:,:,4:6)=1i*e_zx(:,:,:,4:6);
e_yz(:,:,:,4:6)=1i*e_yz(:,:,:,4:6);

[gama,exy,ezx,eyz]=calmode(n,z,neff,pol,k0,z_xy,z_zx,z_yz);
exy(:,:,1)=-exy(:,:,1);ezx(:,:,1)=-ezx(:,:,1);eyz(:,:,1)=-eyz(:,:,1);% orientation

if ~isempty(L);% bessel
%LLL=L(1):L(2);
switch size(L,2);% modif 4 2011 puis 4 2015
case 1;LLL=L;
case 2;LLL=L(1):L(2);
case 3;LLL=L(1):L(2):L(3);
otherwise;LLL=retelimine(L); 
end;
LLL=retcolonne(LLL,1);
f=zeros(1,length(LLL));
ee_xy=enroule(exy,x_xy,y_xy,-LLL,gama,pol);
ee_zx=enroule(ezx,x_zx,y_zx,-LLL,gama,pol);
ee_yz=enroule(eyz,x_yz,y_yz,-LLL,gama,pol);
for iL=1:length(LLL);L=LLL(iL);
f(iL)=-.25i*sqrt(gama/(2*pi))*exp(.5i*(L-.5)*pi)*(diff(ret_int_lorentz(e_xy,ee_xy{iL},wx_xy,wy_xy,[]))+diff(ret_int_lorentz(e_zx,ee_zx{iL},wx_zx,[],wz_zx))+diff(ret_int_lorentz(e_yz,ee_yz{iL},[],wy_yz,wz_yz)));
end;
amp=f*exp(i*LLL.'*retcolonne(teta,1));

else;% cartesien
LLL=[];
amp=zeros(size(teta));
for ii=1:length(teta);
ee_xy=expand_2D(exy,x_xy,y_xy,z_xy,pi+teta(ii),gama,pol);
ee_zx=expand_2D(ezx,x_zx,y_zx,z_zx,pi+teta(ii),gama,pol);
ee_yz=expand_2D(eyz,x_yz,y_yz,z_yz,pi+teta(ii),gama,pol);
amp(ii)=(diff(ret_int_lorentz(e_xy,ee_xy,wx_xy,wy_xy,[]))+diff(ret_int_lorentz(e_zx,ee_zx,wx_zx,[],wz_zx))+diff(ret_int_lorentz(e_yz,ee_yz,[],wy_yz,wz_yz)));
end;
amp=-sqrt(gama/(32*pi))*exp(.25i*pi)*amp;
f=[];
end;% Bessel ou cartesien
% modes à l'origine et declonage
exy0=expand_2D(exy,0,0,z_xy,0,gama,pol);exy0(:,:,:,4:6)=-1i*exy0(:,:,:,4:6);
ezx0=expand_2D(ezx,0,0,z_zx,0,gama,pol);ezx0(:,:,:,4:6)=-1i*ezx0(:,:,:,4:6);
eyz0=expand_2D(eyz,0,0,z_yz,0,gama,pol);eyz0(:,:,:,4:6)=-1i*eyz0(:,:,:,4:6);

angles=struct('f',f,'amp',amp,'F',abs(amp).^2,'L',LLL,'exy',exy0,'ezx',ezx0,'eyz',eyz0,'neff',gama/k0,'Flux_poynting',Flux_poynting);
if isfield(parm,'test') & parm.test==1;
if isfield(parm,'Scale');Scale=parm.Scale;else;Scale=1;end;
trace_boite_2D(  e_xy,x_xy*Scale,y_xy*Scale,z_xy*Scale,    e_zx,x_zx*Scale,y_zx*Scale,z_zx*Scale,    e_yz,x_yz*Scale,y_yz*Scale,z_yz*Scale ,'');		
end;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function e_cartesien=expand_1D(e,x,gama,sens);
e_cartesien=zeros(size(e,1),length(x),3);
if sens==1;
prv=retcolonne(exp(1i*gama*x),1);
else;
prv=retcolonne(exp(-1i*gama*x),1);e(:,:,3)=-e(:,:,3);
end;
for ii=1:3;e_cartesien(:,:,ii)=e(:,:,ii)*prv;end;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function e_cartesien=expand_2D(e,x,y,z,teta,gama,pol);
e_cartesien=zeros(length(z),length(x)*length(y),6);
cost=cos(teta);sint=sin(teta);
[X,Y]=ndgrid(x,y);X=X(:);Y=Y(:);
prv=exp(1i*gama*(cost*X+sint*Y)).';
if pol==0;
e_cartesien(:,:,1)=(-e(:,:,1)*sint)*prv;
e_cartesien(:,:,2)=(e(:,:,1)*cost)*prv;
e_cartesien(:,:,4)=(e(:,:,2)*cost)*prv;
e_cartesien(:,:,5)=(e(:,:,2)*sint)*prv;
e_cartesien(:,:,6)=e(:,:,3)*prv;
else;
e_cartesien(:,:,1)=(e(:,:,2)*cost)*prv;
e_cartesien(:,:,2)=(e(:,:,2)*sint)*prv;
e_cartesien(:,:,3)=e(:,:,3)*prv;
e_cartesien(:,:,4)=(-e(:,:,1)*sint)*prv;
e_cartesien(:,:,5)=(e(:,:,1)*cost)*prv;
end;
e_cartesien=retreshape(e_cartesien,length(z),length(x),length(y),6);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function e_bessel=enroule(e,x,y,L,gama,pol);
e_bessel=cell(1,length(L));
[e_bessel{:}]=deal(zeros(size(e,1),length(x)*length(y),6));
[X,Y]=ndgrid(x,y);X=X(:);Y=Y(:);
R=sqrt(X.^2+Y.^2);
Teta=atan2(Y,X);cosTeta=retdiag(cos(Teta));sinTeta=retdiag(sin(Teta));
[LL,k,kL]=retelimine([L-1,L,L+1]);
[RR,k,kR]=retelimine(R);
J=retbessel('j',LL,gama*RR);
%J=.5*retbessel('h',LL,2,gama*RR);% test
kL=retreshape(kL,[],3);
clear R X Y
for iL=1:length(L);
if pol==0;
e_bessel{iL}(:,:,1)=-.5*e(:,:,1)*(J(kR,kL(iL,3))+J(kR,kL(iL,1))).';
e_bessel{iL}(:,:,2)=.5i*e(:,:,1)*(J(kR,kL(iL,3))-J(kR,kL(iL,1))).';
e_bessel{iL}(:,:,4)=.5i*e(:,:,2)*(J(kR,kL(iL,3))-J(kR,kL(iL,1))).';
e_bessel{iL}(:,:,5)=.5*e(:,:,2)*(J(kR,kL(iL,3))+J(kR,kL(iL,1))).';
e_bessel{iL}(:,:,6)=e(:,:,3)*J(kR,kL(iL,2)).';
else;
e_bessel{iL}(:,:,1)=.5i*e(:,:,2)*(J(kR,kL(iL,3))-J(kR,kL(iL,1))).';
e_bessel{iL}(:,:,2)=.5*e(:,:,2)*(J(kR,kL(iL,3))+J(kR,kL(iL,1))).';
e_bessel{iL}(:,:,3)=e(:,:,3)*J(kR,kL(iL,2)).';
e_bessel{iL}(:,:,4)=-.5*e(:,:,1)*(J(kR,kL(iL,3))+J(kR,kL(iL,1))).';
e_bessel{iL}(:,:,5)=.5i*e(:,:,1)*(J(kR,kL(iL,3))-J(kR,kL(iL,1))).';
end;
for ii=1:6;e_bessel{iL}(:,:,ii)=e_bessel{iL}(:,:,ii)*retdiag(exp(1i*L(iL)*Teta));end;
[e_bessel{iL}(:,:,1),e_bessel{iL}(:,:,2)]=deal(e_bessel{iL}(:,:,1)*cosTeta-e_bessel{iL}(:,:,2)*sinTeta,e_bessel{iL}(:,:,1)*sinTeta+e_bessel{iL}(:,:,2)*cosTeta);% Ex Ey
[e_bessel{iL}(:,:,4),e_bessel{iL}(:,:,5)]=deal(e_bessel{iL}(:,:,4)*cosTeta-e_bessel{iL}(:,:,5)*sinTeta,e_bessel{iL}(:,:,4)*sinTeta+e_bessel{iL}(:,:,5)*cosTeta);% Hx Hy

e_bessel{iL}=retreshape(e_bessel{iL},size(e,1),length(x),length(y),6);
end; % iL
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [gama,e1,e2,e3]=calmode(n,y,neff,pol,k0,y1,y2,y3);
% calcul du mode 0D et du champ normalisé pour divers cotes
y=y(:);n=n(:);hrab_min=.02*pi/real(k0);% attention hrab ne doit pas etre nul pour la normalisation du mode
hrabh=hrab_min;hrabb=hrab_min;
hrabh=max(hrabh,(y1(2)-y(1))*1.1);hrabb=max(hrabb,(y(end)-y1(1))*1.1);
if nargin>6;hrabh=max(hrabh,(y2(2)-y(1))*1.1);hrabb=max(hrabb,(y(end)-y2(1))*1.1);end;
if nargin>7;hrabh=max(hrabh,(y3(2)-y(1))*1.1);hrabb=max(hrabb,(y(end)-y3(1))*1.1);end;
[gama,kmax,vm,er,emode,o,ymode]=retmode(pol,n,-diff(y),k0*neff,[],[hrabh,hrabb,0],[],k0);
tab=[[hrabh;-diff(y);hrabb],n];
init={pol,gama,k0};
sh=retb(init,n(1),1);sb=retb(init,n(end),-1);
% modif 12 2014
e1=retchamp(init,tab,sh,sb,[1,0],{ymode{1}(2:end-1)});
e11=retchamp(init,tab,sh,sb,[0,1],{ymode{1}(2:end-1)});
fac=retcolonne(e1(:,:,1:2))\retcolonne(emode{1}(2:end-1,:,1:2));
ffac=retcolonne(e11(:,:,1:2))\retcolonne(emode{1}(2:end-1,:,1:2));
if retcompare(e1(:,:,1:2)*fac,emode{1}(2:end-1,:,1:2)) < retcompare(e11(:,:,1:2)*ffac,emode{1}(2:end-1,:,1:2))
fac=[fac,0];else;fac=[0,ffac];end;
e1=retchamp(init,tab,sh,sb,fac,{y1-y(end)+hrabb});e1(:,:,2:3)=i*e1(:,:,2:3);% clonage
if nargin>6;e2=retchamp(init,tab,sh,sb,fac,{y2-y(end)+hrabb});e2(:,:,2:3)=i*e2(:,:,2:3);end;
if nargin>7;e3=retchamp(init,tab,sh,sb,fac,{y3-y(end)+hrabb});e3(:,:,2:3)=i*e3(:,:,2:3);end;
% 
% 	e1=retchamp(init,tab,sh,sb,[1,0],{ymode{1}(2:end-1)});
% 	% e1=retchamp(init,tab,sh,sb,[0,1],{ymode{1}(2:end-1)});
% 	fac=retcolonne(e1(:,:,1:2))\retcolonne(emode{1}(2:end-1,:,1:2));
% 	% retcompare(e1(:,:,1:2)*fac,emode{1}(2:end-1,:,1:2)),% pour normalisation
% 	% figure;plot(ymode{1}(2:end-1),real(squeeze(fac*e1(:,1,:))),'-',ymode{1},real(squeeze(emode{1}(:,1,:))),'--',ymode{1}(2:end-1),imag(squeeze(fac*e1(:,1,:))),'-',ymode{1},imag(squeeze(emode{1}(:,1,:))),'--')
% 	e1=retchamp(init,tab,sh,sb,[fac,0],{y1-y(end)+hrabb});e1(:,:,2:3)=i*e1(:,:,2:3);% clonage
% 	if nargin>6;e2=retchamp(init,tab,sh,sb,[fac,0],{y2-y(end)+hrabb});e2(:,:,2:3)=i*e2(:,:,2:3);end;
% 	if nargin>7;e3=retchamp(init,tab,sh,sb,[fac,0],{y3-y(end)+hrabb});e3(:,:,2:3)=i*e3(:,:,2:3);end;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function angles=caldop_modes_popov(n,z,L,teta, e_x,x_x,z_x,wx_x,    e_z,x_z,z_z,wz_z, k0,sym,parm);
[neff,pol,n]=deal(n(end-1),n(end),n(1:end-2));
if nargin<15;parm=struct('');end
% clonage des champs
e_x(:,:,:,4:6)=1i*e_x(:,:,:,4:6);
e_z(:,:,:,4:6)=1i*e_z(:,:,:,4:6);
[gama,ex,ez]=calmode(n,z,neff,pol,k0,z_x,z_z);
ex(:,:,1)=-ex(:,:,1);ez(:,:,1)=-ez(:,:,1);% orientation
ee_x=enroule(ex,x_x,0,-L,gama,pol);
ee_z=enroule(ez,x_z,0,-L,gama,pol);
f=-.25i*sqrt(gama/(2*pi))*exp(.5i*(L-.5)*pi)*(diff(ret_int_lorentz(e_x,ee_x{1},(2*pi*x_x).*wx_x,1,[]))+ret_int_lorentz(e_z,ee_z{1},[],2*pi*x_z,wz_z));
teta=teta(:);
switch sym
case 1;amp=f*cos(L*teta);
%case -1;amp=f*sin(L*teta);
case -1;amp=1i*f*sin(L*teta);% modif 1 2015
otherwise;amp=f*exp(1i*L*teta);
end;

ex(:,:,2:3)=-1i*ex(:,:,2:3);ez(:,:,2:3)=-1i*ez(:,:,2:3);% declonage
angles=struct('teta',teta,'f',f,'amp',amp,'F',abs(amp).^2,'ex',ex,'ez',ez,'neff',gama/k0);
if isfield(parm,'test') & parm.test==1;
trace_boite_popov(   e_x,x_x,z_x,wx_x,   e_z,x_z,z_z,wz_z, '')
end;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%  POPOV  %%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%

function angles=caldop_popov(n,z,u,v,  e_x,x_x,z_x,wx_x,    e_z,x_z,z_z,wz_z, k0,sens,L,sym,parm);
if nargin<16;sym=0;end;
if nargin<17;parm=[];end;

defaultopt=struct('uvmesh',0,'clone',0,'bloch',0);
if isempty(parm);parm=defaultopt;end;
uvmesh=retoptimget(parm,'uvmesh',defaultopt,'fast');
clone=retoptimget(parm,'clone',defaultopt,'fast');
bloch=isfield(parm,'tab');
if ~uvmesh;
nu=length(u);nv=length(v);nuv=nu*nv;
[u,v]=ndgrid(u,v);
end;
if isempty(z);z=0;n=[n,n];end;
u=u(:);v=v(:);
% clonage des champs
if clone==0;
e_x(:,:,:,4:6)=1i*e_x(:,:,:,4:6);
e_z(:,:,:,4:6)=1i*e_z(:,:,:,4:6);
end;
Flux_poynting=diff(retpoynting({e_x},[0,0,1],x_x.*wx_x,2*pi,[]))+retpoynting({e_z},[1,0,0],[],2*pi*x_z,wz_z);

if sens==1;inc=[0,1];n_op=n(1);z_op=z(1);else;inc=[1,0];n_op=n(end);z_op=z(end);end;
khi=sqrt(max(eps,(k0*n_op)^2-u.^2-v.^2));
delt=0;
	if sens==1;
	K=[u,v,khi];K=K./repmat(sqrt(sum((K.^2),2)),1,size(K,2));%vecteurs d'onde normalises
	V=[-K(:,2),K(:,1),zeros(size(K,1),1)];
	ii=find(abs(K(:,1).^2+K(:,2).^2)<100*eps);V(ii,1)=-sin(delt);V(ii,2)=cos(delt);%incidence normale: delt valeur par defaut de delt
	V=V./repmat(sqrt(sum((V.^2),2)),1,size(V,2));
	U=[K(:,3).*V(:,2),-K(:,3).*V(:,1),K(:,2).*V(:,1)-K(:,1).*V(:,2)];%U= V vect.  K
	else;
	K=[u,v,-khi];K=K./repmat(sqrt(sum((K.^2),2)),1,size(K,2));%vecteurs d'onde normalises
	V=[-K(:,2),K(:,1),zeros(size(K,1),1)];
	ii=find(abs(K(:,1).^2+K(:,2).^2)<100*eps);V(ii,1)=-sin(delt);V(ii,2)=cos(delt);%incidence normale: delt valeur par defaut de delt
	V=V./repmat(sqrt(sum((V.^2),2)),1,size(V,2));
	U=[K(:,3).*V(:,2),-K(:,3).*V(:,1),K(:,2).*V(:,1)-K(:,1).*V(:,2)];%U= V vect.  K
	end
delta=atan2(-V(:,1),V(:,2));cos_delta=cos(delta);sin_delta=sin(delta);
[uv_reduit,prv,kkuv]=retelimine(sqrt(u.^2+v.^2),-1);

Z=[z_x(:);z_z(:)];% tous les Z où on veut le champ
[Z_reduit,prv,kkZ]=retelimine(Z,-1);
% calcul des Bessels (le 'r' de l'integration intervient ici)
if sym==0;
JL=retbessel('j',[L-1,L,L+1],retcolonne(uv_reduit*([x_x(:);x_z(1)].')));% on met en dernier le Z_x qui est constant
JLM1=retreshape(JL(:,1),length(uv_reduit),length(x_x)+1)*retdiag([x_x(:);x_z(1)]);
JLP1=retreshape(JL(:,3),length(uv_reduit),length(x_x)+1)*retdiag([x_x(:);x_z(1)]);
JL=retreshape(JL(:,2),length(uv_reduit),length(x_x)+1)*retdiag([x_x(:);x_z(1)]);
else;
JML=retbessel('j',[L-1,L,L+1,-L-1,-L,-L+1],retcolonne(uv_reduit*([x_x(:);x_z(1)].')));% on met en dernier le Z_x qui est constant
JLM1=retreshape(JML(:,1),length(uv_reduit),length(x_x)+1)*retdiag([x_x(:);x_z(1)]);
JLP1=retreshape(JML(:,3),length(uv_reduit),length(x_x)+1)*retdiag([x_x(:);x_z(1)]);
JL=retreshape(JML(:,2),length(uv_reduit),length(x_x)+1)*retdiag([x_x(:);x_z(1)]);
JMLM1=retreshape(JML(:,4),length(uv_reduit),length(x_x)+1)*retdiag([x_x(:);x_z(1)]);
JMLP1=retreshape(JML(:,6),length(uv_reduit),length(x_x)+1)*retdiag([x_x(:);x_z(1)]);
JML=retreshape(JML(:,5),length(uv_reduit),length(x_x)+1)*retdiag([x_x(:);x_z(1)]);
end;
if isempty(z);z=mean(Z_reduit);n=[n(1),n(1)];end;
hh=max(0,max(Z_reduit)-max(z))+1.e-6;hb=max(0,min(z)-min(Z_reduit))+1.e-6;


tab=[[hh;-diff(z(:));hb],n(:)];
tab0=tab;tab0(:,2)=n_op;

e=cell(1,2);
for kpol=1:2;
init={2*kpol-2,-uv_reduit,k0};
if bloch & sens==-1;[sh,prv]=retbloch(init,retcouche(init,parm.tab),sum(parm.tab(:,1)));else;sh=retb(init,n(1),1);end;
if bloch & sens==1; [prv,sb]=retbloch(init,retcouche(init,parm.tab),sum(parm.tab(:,1)));else;sb=retb(init,n(end),-1);end;
sh_inc=retb(init,n_op,1);sb_inc=retb(init,n_op,-1);
e_inc=retchamp(init,tab0,sh_inc,sb_inc,inc,{z_op-z(end)+hb});e_inc(:,:,2:3,:,:)=1i*e_inc(:,:,2:3,:,:);% clonage
e{kpol}=retchamp(init,tab,sh,sb,inc,{Z_reduit-z(end)+hb});e{kpol}(:,:,2:3,:,:)=1i*e{kpol}(:,:,2:3,:,:);% clonage


for ii=1:size(e{kpol},5);e{kpol}(:,:,:,:,ii)=e{kpol}(:,:,:,:,ii)/e_inc(:,:,1,:,ii);end;% normalisation E=1
end;  % kpol

EE=zeros(length(u),2);
HH=zeros(length(u),2);

for kpol=1:2;
for ii=1:length(u);
eee_x=cale_popov(x_x,e{kpol}(kkZ(1:length(z_x)),:,:,:,kkuv(ii)),size(e_x),cos_delta(ii),sin_delta(ii),kpol,L,delta(ii),JLM1(kkuv(ii),1:end-1),JL(kkuv(ii),1:end-1),JLP1(kkuv(ii),1:end-1));
eee_z=cale_popov(x_z,e{kpol}(kkZ(1+length(z_x):length(z_x)+length(z_z)),:,:,:,kkuv(ii)),size(e_z),cos_delta(ii),sin_delta(ii),kpol,L,delta(ii),JLM1(kkuv(ii),end),JL(kkuv(ii),end),JLP1(kkuv(ii),end));
if sym~=0;
e_xx=e_x;e_zz=e_z;e_xx(:,:,:,[2,4,6])=-e_xx(:,:,:,[2,4,6]);e_zz(:,:,:,[2,4,6])=-e_zz(:,:,:,[2,4,6]);
eee_xx=cale_popov(x_x,e{kpol}(kkZ(1:length(z_x)),:,:,:,kkuv(ii)),size(e_x),cos_delta(ii),sin_delta(ii),kpol,-L,delta(ii),JMLM1(kkuv(ii),1:end-1),JML(kkuv(ii),1:end-1),JMLP1(kkuv(ii),1:end-1));
eee_zz=cale_popov(x_z,e{kpol}(kkZ(1+length(z_x):length(z_x)+length(z_z)),:,:,:,kkuv(ii)),size(e_z),cos_delta(ii),sin_delta(ii),kpol,-L,delta(ii),JMLM1(kkuv(ii),end),JML(kkuv(ii),end),(-1)^(L+1)*JMLP1(kkuv(ii),end));
end;
switch sym;
case 0;EH=diff(ret_int_lorentz(e_x,eee_x,wx_x,1,[]))+ret_int_lorentz(e_z,eee_z,[],1,wz_z);
case 1;EH=.5*(diff(ret_int_lorentz(e_x,eee_x,wx_x,1,[]))+ret_int_lorentz(e_z,eee_z,[],1,wz_z)+diff(ret_int_lorentz(e_xx,eee_xx,wx_x,1,[]))+ret_int_lorentz(e_zz,eee_zz,[],1,wz_z));
case -1;EH=.5*(diff(ret_int_lorentz(e_x,eee_x,wx_x,1,[]))+ret_int_lorentz(e_z,eee_z,[],1,wz_z)-diff(ret_int_lorentz(e_xx,eee_xx,wx_x,1,[]))-ret_int_lorentz(e_zz,eee_zz,[],1,wz_z)); 
end;

if kpol==1;EE(ii,2)=k0/(8i*pi^2*khi(ii))*EH;else;HH(ii,2)=k0*n_op^2/(8i*pi^2*khi(ii))*EH;end;
end;%ii
end;%kpol

mu=k0;ep=k0*n_op^2;
HH(:,1)=-1i*(k0*n_op/mu)*EE(:,2);% -1 à cause de l'orientation
EE(:,1)=-1i*(k0*n_op/ep)*HH(:,2);

HH=-1i*HH;% declonage
% analyse des polarisations
psi=.5*atan2(real(EE(:,1)).*real(EE(:,2))+imag(EE(:,1)).*imag(EE(:,2)),.5*(abs(EE(:,1)).^2-abs(EE(:,2)).^2));
EEE(:,1)=EE(:,1).*cos(psi)+EE(:,2).*sin(psi);EEE(:,2)=-EE(:,1).*sin(psi)+EE(:,2).*cos(psi);
ii=find(abs(EEE(:,2))>abs(EEE(:,1)));psi(ii)=psi(ii)+pi/2;[EEE(ii,1),EEE(ii,2)]=deal(EEE(ii,2),-EEE(ii,1));
ii=find(psi<=-pi/2+100*eps);psi(ii)=psi(ii)+pi;EEE(ii,:)=-EEE(ii,:);
HHH(:,1)=HH(:,1).*cos(psi)+HH(:,2).*sin(psi);HHH(:,2)=-HH(:,1).*sin(psi(:))+HH(:,2).*cos(psi);
ii=find(abs(psi-pi/2)<100*eps);psi(ii)=pi/2;ii=find(abs(psi)<100*eps);psi(ii)=0;
teta=acos(K(:,3));

F=2*pi^2*(k0*n_op)^2*retdiag(cos(teta).^2)*real(EE(:,1).*conj(HH(:,2))-EE(:,2).*conj(HH(:,1)));
e=full([retdiag(EE(:,1))*U+retdiag(EE(:,2))*V,retdiag(HH(:,1))*U+retdiag(HH(:,2))*V]);
if sens>0;theta=teta;else;theta=pi-teta;end;
E=-2i*pi*n_op*k0*retdiag(cos(theta))* [ (retdiag(EE(:,1))*U+retdiag(EE(:,2))*V),(retdiag(HH(:,1))*U+retdiag(HH(:,2))*V)];

if ~uvmesh;
e=retreshape(e,nu,nv,[]);
EE=retreshape(EE,nu,nv,[]);
HH=retreshape(HH,nu,nv,[]);
EEE=retreshape(EEE,nu,nv,[]);
HHH=retreshape(HHH,nu,nv,[]);
F=retreshape(full(F),nu,nv,[]);
end;
angles=struct('teta',teta,'theta',theta,'delta',delta,'psi',psi,'k',K,'u',U,'v',V,'e',-E,'E',E,'EE',EE,'HH',HH,'EEE',EEE,'HHH',HHH,'F',full(F),'origine',z_op,'Flux_poynting',Flux_poynting);
if isfield(parm,'test') & parm.test==1;
trace_boite_popov( e_x,x_x,z_x,wx_x,   e_z,x_z,z_z,wz_z ,'');
end;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function eee=cale_popov(X,e,sz,cos_delta,sin_delta,kpol,L,delta,JLM1,JL,JLP1);
delta=delta+pi;
eee=zeros(sz);
if kpol==1;
eee(:,:,:,2)=-retreshape(repmat(e(:,1,1),1,length(X)),sz(1:3));
eee(:,:,:,4)=retreshape(repmat(e(:,1,2),1,length(X)),sz(1:3));
eee(:,:,:,6)=retreshape(repmat(e(:,1,3),1,length(X)),sz(1:3));
else;
eee(:,:,:,5)=-retreshape(repmat(e(:,1,1),1,length(X)),sz(1:3));
eee(:,:,:,1)=retreshape(repmat(e(:,1,2),1,length(X)),sz(1:3));
eee(:,:,:,3)=retreshape(repmat(e(:,1,3),1,length(X)),sz(1:3));
end;

[eee(:,:,:,1),eee(:,:,:,2)]=deal(eee(:,:,:,1)*cos_delta-eee(:,:,:,2)*sin_delta,eee(:,:,:,1)*sin_delta+eee(:,:,:,2)*cos_delta);
[eee(:,:,:,4),eee(:,:,:,5)]=deal(eee(:,:,:,4)*cos_delta-eee(:,:,:,5)*sin_delta,eee(:,:,:,4)*sin_delta+eee(:,:,:,5)*cos_delta);

JLP1=(pi*exp(1i*((L+1)*delta+L*pi/2)))*JLP1;
JL=(pi*exp(1i*L*(delta+pi/2)))*JL;
JLM1=(pi*exp(1i*((L-1)*delta+L*pi/2)))*JLM1;
if size(eee,1)==2;JL=[JL(:).';JL(:).'];JLP1=[JLP1(:).';JLP1(:).'];JLM1=[JLM1(:).';JLM1(:).'];% pour e_x
else;JL=repmat(JL,1,size(eee,1));JLP1=repmat(JLP1,1,size(eee,1));JLM1=repmat(JLM1,1,size(eee,1)); % pour e_z
end;
JL=retreshape(JL,size(eee(:,:,:,1)));JLP1=retreshape(JLP1,size(eee(:,:,:,1)));JLM1=retreshape(JLM1,size(eee(:,:,:,1)));
[eee(:,:,:,1),eee(:,:,:,2),eee(:,:,:,3)]=deal(...
(1i*eee(:,:,:,1)+eee(:,:,:,2)).*JLP1-(1i*eee(:,:,:,1)-eee(:,:,:,2)).*JLM1,...
(-eee(:,:,:,1)+1i*eee(:,:,:,2)).*JLP1-(eee(:,:,:,1)+1i*eee(:,:,:,2)).*JLM1,...
2*eee(:,:,:,3).*JL);
[eee(:,:,:,4),eee(:,:,:,5),eee(:,:,:,6)]=deal(...
(1i*eee(:,:,:,4)+eee(:,:,:,5)).*JLP1-(1i*eee(:,:,:,4)-eee(:,:,:,5)).*JLM1,...
(-eee(:,:,:,4)+1i*eee(:,:,:,5)).*JLP1-(eee(:,:,:,4)+1i*eee(:,:,:,5)).*JLM1,...
2*eee(:,:,:,6).*JL);

