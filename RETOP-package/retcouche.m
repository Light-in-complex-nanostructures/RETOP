function [a,d,w,maillage]=retcouche(init,w,c,li,maillage);
% function [a,d,w]=retcouche(init,w,c,li);
%
%  calcul du 'descripteur de texture' a, utilisee par retb, retc, retchamp pour calculer les champs
%  si on n'utilise que retb il suffit de prendre c=-1 (ce qui diminue la memoire et le calcul) 
%  si on veut tracer les champs c=1 mais plus long et gourmand en memoire
%  si on ne veut pas tracer les champs c=0 (option par defaut) 
%  init calcule par retinit
%  w est le maillage de texture decrivant l'objet calculee par retu
%  c,li parametres
%
% si c=-inf matrices vides pour tests
%
% POUR LES DIELECTRIQUES
%........................
%
%  pour calculer 'proprement' les champs prendre c=1,mais ceci alourdit les calculs 
%     ( pour les sources  en 2 D si on veut calculer l'energie diffractee par la source il faut prendre c=1
%         ou bien c=[1,x,y] ou x y sont les coordonnees de la source      )
%
%  d:valeurs propres
%  la  matrice m est sous la forme  m=[0,a1;a2,0]=.5*[q,q;p*d,-p*d]*[d,0;0,-d]*[qq,(1/d)*pp;qq,-(1/d)*pp]  
%  en 1D a={    p,pp,q,qq,              d,                      a3,a1,         w                 parm};
%           base vecteurs propres  valeurs propres     matrices servant au calcul des champs   
%                                                        (vides si c=-1ou 0)
%
%  en 2D a={    p,pp,q,qq,               d,              ez,hz,fex,fhx,fey,fhy,  w            pas,       parm};
%           base vecteurs propres  valeurs propres     matrices servant au calcul des champs   pas      
%                                                        (vides si c=-1ou 0)
%   parm=struct('dim',dim,'genre',2,'type',1,'sog',sog);
%   type:  1 dielectriques  ,2 metaux infiniment conducteurs ,3 modes de bloch non symetriques,4 elements finis,5 textures inclinees
%   genre 2  descripteur de texture
%
%  ce programme necessite en general la DIAGONALISATION d'une matrice n*n 
%  (sauf pour les metaux et les milieux homogenes isotropes)
% 
%  en 2D 
%
%  si on prend imag(c)~=0 ,les matrices sont stockées sur disque(ne pas oublier de faire retio a la fin pour liberer les fichiers)
%            si de plus imag(c)=-2  en 2 D  on utilise une bibliotheque
%            pour effacer la bibliotheque :retcouche;  
%            pour recuperer cette bibliotheque :bibli=retcouche;  
%            pour regenerer cette bibliotheque :retcouche(bibli);
%                      attention si bibli contient des fichiers a ce qu'ils ne soient pas effaces  
%
% si w est mis en sortie il devient: w={w{1},w{2},w{3},eepz,mmuz,muy,mux,epy,epx,hx,ez,hz,fex,fhx,fey,fhy}
% ce qui permet de stocker les divers tableaux eepz,...
% dans un  nouvel appel ils ne sont plus recalcules
%  si de plus c=-2 on ne fait pas la diagonalisation ( utile a retbrillouin qui n'utilise que eepz,mmuz,muy,mux,epy,epx)  
% ceci permet d'eviter des calculs inutiles pour des objets compliques quand lambda varie 
%  ATTENTION les eps ne doivent pas avoir varié entre temps...
% si de plus imac(c)~=0  eepz,...etc sont stockés sur disque  (ne pas oublier de faire retio a la fin )%
%
%   li(12):parametres(facultatif) pour la methode de li (1 -1 0 ou 2) par defaut li=zeros(1,12) ou 12 fois li(1) si 1 seul terme. 
%   si 1  X->Y puis Y->X 
%   si -1 Y->X puis X->Y 
%   si 2 pas methode de Li ( milieux continus)
%   si 0 le choix du sens depend de la continuitee des fonctions
%    li         | 1   |  2  | 3 | 4 | 5 | 6 |            li        | 7 | 8 | 9 | 10| 11| 12|
%               -----------------------------                      -------------------------
%    matrices   |eepz | mmuz|muy|mux|epy|epx|   calcul des champs: |epz|muz|fex|fhx|fey|fhy|
% en conique, si les textures sont invariantes en y il est preferable de prendre li=1
% 
% POUR LES TEXTURES INCLINEES
% le maillage de texture w décrivant l'objet, calculée par retu doit être remplacé par:
% en 1D   {w,struct('teta',teta,'fac',fac,'Maystre',Maystre)}
% en 2D   {w,struct('teta',teta,'phi',phi,'sens',sens,'fac',fac,'Maystre',Maystre)}
%   teta, phi,angles en radians,
%   Maystre angle de la coupure en radians (facultatif sinon coupure fonction de teta)
%   sens = 1 ou -1 suivant que l'on travaille en x ou en y
%   fac=[1+.01i ,1+.1i] pour tri des modes par exemple on diagonalise pour k0 k0*(1+.01i) et on extrapole à k0*(1+.1i)
%   sens et fac sont facultatifs (par défaut sens = 1 fac=[])
%
%
% 1D
% function [a,dd]=retcouche(init,X,ep,teta);
%  milieu periodique invariant 
% dans une direction inclinee d'un angle teta(degres) sur y
% ep valeurs de [epsilon;mux;1/muy]
% X={x,wx} points et poids de Gauss 
% maillage de texture: u={{x,wx},ep,teta,struct('type',5)}
% Attention si l'inclinaison est non nulle pas de changement de coordonnees !
% EXEMPLE reseau sinusoidal 
%   mm=40;init=retinit(d,[-mm,mm],beta,sym);
%   dnn=1.e-2;inclinaison=10;
% 	ntr=max(2*mm+1,5);[x,w]=retgauss(0,d,ntr,-4);nn=nn0+dnn*cos(2*pi*x/d);% echantillonage sinusoide
% 	ep=retep(nn,pol,k0);x={x,w};% ep valeurs de [epsilon;mux;1/muy]
% 	u={x,ep,inclinaison*pi/180,struct('type',5)};
% 
% 2D
% function [a,dd]=retcouche(init,X,Y,ep,teta,phi);
% X={x,wx} ,Y={y,wy} points et poids de Gauss 
% maillage de texture: u={{x,wx},{y,wy},ep,teta,phi,struct('type',5)}
% 
% 
% RECAPITULATIF pour le parametre c
% |   c   |retchamp|retc|retb|rets|retbrillouin|retneff|cpu                             |imag(c)|
% ----------------------------------------------------------                            ---------
% |   1   | oui    |oui |oui |oui |oui         |oui    |1      calcul complet           | ~=0   | stockage sur disque
% |[1,x,y]|approché|oui |oui |oui |oui         |oui    |1      champ sur la source      |    -2 | bibliotheque(obsolete)  
% |   0   |approché|oui |oui |non |oui         |oui    |.99    le champ est calculé approximativement
% |  -1   |non     |non |oui |non |oui         |oui    |.9     seulement conditions aux limites(retb) 
% |  -2   |non     |non |non |non |oui         |non    |.01    pas diagonalisation 
% |  -3   |non     |non |non |non |non         |oui    |.5     seulement les valeurs propres  
% | -inf  |non     |non |non |non |non         |non    |0      matrices vides pour tests  
%
%..............................................
% POUR LES METAUX INFINIMENT CONDUCTEURS
%..............................................
%  function [a,d]=retcouche(init,w,c);
%  calcul du descripteur de texture a utilisee par  retc pour le metal infiniment conducteur
%  a={sparse(g1),sparse(g2),dd,m,n,champ,www};
%  dd:valeurs propres
%  c:0 ou 1 si on veut calculer les champs
%
% 
%   ELEMENTS FINIS
%..............................................
%  [a,m,me,maillage]=retcouche(init,w,parm,k0,maillage);
%  w est un maillage de tronçon obtenu par retautomatique
%  k0=2*pi/ld normalisation des eps mu si les unites sont metriques (1 par defaut)
%    (les eps de maillage.eps sont calcules pour k0=1 et pour chaque calcul multiplies par k0)
%   si maillage existe en entree et n'est pas vide,le maillage n'est pas refait (gain de temps si par exemple seul,ld varie) 
%     (par exemple: maillage=[];for ld= ...;[a,m,me,maillage]=retcouche(init,w,parm,2*pi/ld,maillage);end; )
%  en sortie:
%   a  descripteur du tronçon  type elements finis
%  m, me   nb de degres de liberte , nb d' elements
%
% parm: parametres du maillage sous forme d'une structure à champs
% parm=struct('h0',.. ,'poids',[] ,'choix', 0,'c',0,'serre',1.3,'tracer',0,.. ,'Pfix',[],'pol',0,'inclusions',{ inclusions },'pasmax',3 )
% h0: pas 'approximatif' du maillage  
% poids :  fonction de variable d'espace x(:,dim) permettant de moduler le maillage. facultatif (par defaut [] )
% 	  le pas est proportionnel à poids(x) avec un facteur de normalisation tel que le nombre d'elements 
% 	  reste le même (surface/h0^2 ou volume/h0^3)
% 
% choix:  type de maillage facultatif ( 0 par defaut: Delaunay )
% 	  en 1D choix=1 : ancien maillage ,choix=0: Delaunay 
% 	  en 2D choix=1 : 5 tetraedres par cube ,choix= 2:  6 tetraedres par cube ,choix=0: Delaunay 
% 
% c=1 preparation du calcul du champ  (c=0 par defaut)
% 	  si la partie imaginaire de c n'est pas nulle,mise d'elements sur fichiers (ne pas oublier de faire retio ..)
%
% serre parametre qui permet d'empecher les tetraedres de traverser les surfaces maillees.1.3 par defaut
% 	  si serre est trop grand le maillage est irregulier au voisinage des arêtes
% 
% pasmax:iterations de Pearson ( par defaut:3 en 1 D ,5 en 2 D)
% 	tracer:
%    en 1D: 0 pas trace, 1 trace du maillage
%    en 2D: 0 pas trace,1:une seule vue,2:une figure par domaine,3:une figure avec tous les domaines separes en subplots,
%          4:vue de l'objet coupé en 8 sur une figure,5:une figure par domaine coupe en 8
%    il est possible de demander plusieurs traces (tracer est un vecteur)
% 	pol  en 2 D, pol  0  elements finis en E  , 2 elements finis en H
%  'Pfix':[x(:),y(:)] ou [x(:),y(:),z(:)] points fixes 
% inclusions: cell array de structures a champs qui permet d' ajouter au maillage des elements nouveaux.
% 	Ces elements 'ecrasent' l'objet defini par w et s'ecrasent mutuellement (le dernier etant prioritaire)
% 	inclusions{ii}=struct('type','polygone','sommets',[xx.',yy.'],'pointindice',retep(nm,pol),'r',[],'centre',[]);
% 	si une seule inclusion inutile de mettre en cell array
% 	ATTENTION: si on definit parm avec l'ordre struct,le cell array doit etre mis lui même entre {  }  
% 
% 		en 1 D:   les champs de inclusions  sont:
% 		'type' :  'cercle'  (en fait ellipse)   ou  'polygone' defini par des points  
% 		'indice',[n;pol], [xi,yi]  ou retep(n,pol) pour definir l'indice de l'inclusion,par l'indice et pol (vecteur colonne)
% 		par ses parametres electromagnetiques ou comme le même indice qu'un point de l'objet
% 		si type= 'cercle',les autres champs sont:
% 			'centre', [x0,y0],
% 			'r' ,   [rx,ry]  (dans le cas d'un cercle une seule valeur suffit) 
% 		si type= 'polygone',l' autre champs est:
% 			'sommets', [x(:),y(:)],  le polygone n'est pas fermé dans le maillage ,
% 			mais l'indice est mis a l'interieur du polygone.Il y a un noeud du maillage a chaque sommet
% 			du polygone, sauf si on ajoute i a ses coordonnees 
% 			(cas de courbes decrites par beaucoup de points)  
% 
% 		en 2 D:   les champs de inclusions  sont:
% 		'direction':[teta,phi] (en radiants)
% 		'centre': [x0,y0,z0]
% 		l'inclusion est definie dans le repere Oxyz,puis subit une rotation autour du centre ,d'angles direction (voir retdeplacement)
%       'Pfix':[x(:),y(:),z(:)] points fixes qui sont transformes par le deplacement ci dessus 
% 		'volume' :0 on ne maille la surface sans definir un domaine avec un indice different  
% 		'type' :  'sphere' (en fait ellipsoide),
% 		'cylindre'(en fait cylindre de bases elliptiques differentes)
% 		'polyedre' (defini comme l'enveloppe convexe d'un ensemble de sommets)  
% 		'indice', n ,[xi,yi,zi]  ou ret2ep(n) pour definir l'indice de l'inclusion par un indice,
% 		ses parametres electromagnetiques ou comme le même indice qu'un point de l'objet
% 		si type= 'sphere',les autres champs sont:
% 			'centre', [x0,y0,z0],
% 			'r' ,   [rx,ry,rz]  (dans le cas d'une sphere une seule valeur suffit) 
% 		si type= 'cylindre',les autres champs sont:
% 			'h' ,  hauteur (si h=0 le cylindre degenere en une ellipse(couvercle=[1,-1]))
% 			'r' ,   [[r_haut_x,r_haut_y];[r_bas_x,r_bas_y]] ou [r_haut_x;r_bas_x] ,ou, r
% 			le cylindre peut etre dégénére en cone ou cylindre 'pincé' sur une face;
% 			'couvercle' : [ch,cb]  ch,cb valant:1 on maille le couvercle et son contour                                              
% 			0 :on ne maille pas le couvercle mais on maille son contour                                              
% 			-1: on  maille ni le couvercle ni son contour                                              
% 		si type= 'polyedre',l' autre champ est:
% 			'sommets' [x(:),y(:),y(:)]
% %% EXEMPLES 2D
% d=[.6,.6];h=.3;nh=1;nb=1.5;n1=1;n2=3.5;nn=5i;centre=d/2;sym=[1,1,centre];mm=6;ld=.6;k0=2*pi/ld;
% unv=retu([d(1),2*h],{[h/2;3*h/2;2*h],[nb;n1;nh]});w=retautomatique(d,[],unv,2*h);
% unv1=retu([d(1),2*h],{[h/2;3*h/2;2*h],[nh;n2;nh]});w1=retautomatique(d,[],unv1,2*h);
% unv0=retu([d(1),2*h],{n1});w0=retautomatique(d,[],unv0,2*h);
% init=retinit(d,[-mm+i,mm,-mm,mm],[],sym);
% % pyramide
% inclusions=struct('type','polyedre','centre',[centre,h/2],'sommets',[[-d(1)/4,-d(2)/4,0];[d(1)/4,-d(2)/4,0];[-d(1)/4,d(2)/4,0];[d(1)/4,d(2)/4,0];[0,0,h]],'indice',nn,'direction',[pi/4,pi/2]);
% [a,Nd]=retcouche(init,w,struct('h0',sqrt(prod(d))/10,'inclusions',{inclusions},'tracer',3),k0);
%% cone
%inclusions=struct('type','cylindre','centre',[centre,h],'r',[[0,0];[d(1)/3,d(2)/3]],'h',h,'indice',nn);
% [a,Nd]=retcouche(init,w,struct('h0',sqrt(prod(d))/10,'inclusions',{inclusions},'tracer',3),k0);
%% cylindre pose verticalement sur un dioptre 
% inclusions=struct('type','cylindre','centre',[centre,h],'r',d(1)/3,'h',h,'indice',nn);
% [a,Nd]=retcouche(init,w,struct('h0',sqrt(prod(d))/10,'inclusions',{inclusions},'tracer',[3,4]),k0);
%% cylindre d'axe OX
% inclusions=struct('type','cylindre','centre',[centre,h],'r',d(1)/3,'h',d(2),'indice',nn,'direction',[0,0]);
% [a,Nd]=retcouche(init,w0,struct('h0',sqrt(prod(d))/10,'inclusions',{inclusions},'tracer',3),k0);
%% cylindre d'axe OY
% inclusions=struct('type','cylindre','centre',[centre,h],'r',d(2)/3,'h',d(1),'indice',nn,'direction',[pi/2,0]);
% [a,Nd]=retcouche(init,w0,struct('h0',sqrt(prod(d))/10,'inclusions',{inclusions},'tracer',3),k0);
%% cylindre 'pincé'
% inclusions=struct('type','cylindre','centre',[centre,h],'r',[[d(1)/6,0];[d(1)/3,d(2)/3]],'h',h,'indice',nn);
% [a,Nd]=retcouche(init,w,struct('h0',sqrt(prod(d))/10,'inclusions',{inclusions},'tracer',3),k0);
%% cone sur cylindre
% inclusions={struct('type','cylindre','centre',[centre,3*h/4],'r',.6*h,'h',h/2,'couvercle',[-1,1],'indice',nn),...
% struct('type','cylindre','centre',[centre,5*h/4],'r',[[0,0];[.6*h,.6*h]],'h',h/2,'couvercle',[0,0],'indice',nn)};
% [a,Nd]=retcouche(init,w,struct('h0',sqrt(prod(d))/10,'inclusions',{inclusions},'tracer',3),k0);
%% ellipsoide dans milieu homogene
%inclusions=struct('type','sphere','centre',[centre,h],'r',[d/3,.8*h],'indice',nn);
% [a,Nd]=retcouche(init,w0,struct('h0',sqrt(prod(d))/10,'inclusions',{inclusions},'tracer',[3,4]),k0);
%% maille d'un cristal photonique
% inclusions={struct('type','cylindre','centre',[centre,h],'r',d(1)/4,'h',h,'couvercle',[0,0],'indice',n1),...
% struct('type','cylindre','centre',[centre+d/2,h],'r',d(1)/4,'h',h,'couvercle',[0,0],'indice',nh)};
% [a,Nd]=retcouche(init,w1,struct('h0',sqrt(prod(d))/10,'inclusions',{inclusions},'tracer',3,'serre',1.4),k0);
%% sphère en partie enfoncée dans un dioptre. Ici le cylindre de hauteur 0 sert a mailler l'intersection de la sphere avec le dioptre
% inclusions={struct('type','sphere','centre',[centre,h],'r',sqrt(2)*h/2,'indice',nn,'direction',[0,pi/2],'h',[]),...
% struct('type','cylindre','centre',[centre,h/2],'r',h/2,'h',0)};
% [a,Nd]=retcouche(init,w,struct('h0',sqrt(prod(d))/10,'inclusions',{inclusions},'tracer',3),k0);
%
%% EXEMPLE 1D
% pol=0;ld=1.5;d=3;init=retinit(d,[-1+i,1],[],[1,0]);w={{3,retu(d,{1,pol})}};
% inclusions={struct('type','polygone','sommets',[[-1;0;1;0;-1],[1.5;2.5;1.5;.5;1.5]],'indice',[3.5;pol]),...
% struct('type','cercle','centre',[0,1.5],'r',[.5,.5],'indice',[1;pol])};
%  poids=@(x) (1.5-retpoids(x,[0,1.5],[.5,.5],d));
% [a_elf,md,me]=retcouche(init,w,struct('h0',.1,'tracer',1,'inclusions',{inclusions},'poids',poids),2*pi/ld); %  maillage elements finis
%
%  See also:RETB,RETC,RETBLOCH,RETCHAMP,RETMAILLAGE,RETDEPLACEMENT,RETTETRAMESH,RETPOIDS,RETHELP_POPOV


persistent bibli;  % utilisation d'une bibliotheque
if nargin==0;if nargout==0;bibli=retio(bibli,-1);else a=bibli;end;return;end;
if nargin==1;bibli=init;end; % on regenere la bibliotheque

if nargin<3;c=0;end;if length(c)>1;xsource=c(2:end);c=c(1);else xsource=[];end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~isstruct(init{end}); % 0 D rapide
%[prv1,prv2,prv3,prv4,a]=retabeles(init{1},w(:,2:end),w(:,1),init{2},init{3},[],0);
[prv1,prv2,prv3,prv4,a]=retabeles(init{1},w(:,2),w(:,1),init{2},init{3},[],0);
return;end;% fin 0 D rapide
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Cercles Aureliens
if ~iscell(w{1}) & iscell(w{end});
[a,d,w]=retcouche_Aurelien(init,w,c);	
return;end;% fin Cercles Aureliens
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% elements finis
if iscell(w{end});
if nargin<4;li=1;end;
if nargin<5;maillage=[];end;
[a,d,w,maillage]=retelf_couche(init,w,c,li,retio(maillage));
return;end;% fin elements finis
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% fonctions continues ou inclinaison (à laisser apres ELF)
if isstruct(w{end});
[a,d]=retincline(init,c,w);
return;end;% fin fonctions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% cylindres Popov
if init{end}.genre==1;
[a,d]=retcouche_popov(init,w,c);
return;end;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin<4;li=0;end;if isempty(li);li=0;end;if length(li)<12;li=li(1)*ones(1,12);end; % doit etre apres elements finis
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if c==-inf; % matrices vides pour tests
if size(w,2)==5 ;% metaux
a=cell(1,7);a{7}=w;
else;            % dielectriques
if init{end}.dim==2; %2D
a=cell(1,14);a{12}=w;
else;
a=cell(1,9);a{8}=w;
end;
end;
d=[];
return;end;  % fin matrices vides    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
if size(w,2)==5 ;% < +++++++ metaux
[a,d]=retmcouche(init,w,real(c));
else;            % < +++++++  dielectriques
if init{end}.dim==2; % < ******************* 2D 
if imag(c)==-2; % utilisation d'une bibliotheque
if isempty(init{10});iinit=[init{2}(:);init{end}.sym(:)];else;iinit=[init{2}(:);init{end}.sym(:);init{10}{1}(:);init{10}{2}(:)];end;
if length(w)==5;ww=[w{1}(:);w{2}(:);w{3}(:);w{4}(:);w{5}(:)];else;ww=[w{1}(:);w{2}(:);w{3}(:)];end; 
if ~isempty(bibli);
if abs(retcompare(iinit,bibli{1},2))<10*eps;    
for ii=2:length(bibli);
if (abs(retcompare(ww,bibli{ii}{1},2))<10*eps)&(real(c)<=bibli{ii}{2});% valeur trouvee
a=bibli{ii}{3};w=bibli{ii}{4};return;end;
end;
else;bibli=retio(bibli,-1);end;
end;
end; 

if ~isempty(init{end}.L);if nargin>3;c=c+i;end;% cylindriquee_radial
if init{end}.granet;pas=init{end}.d;x=w{1};y=w{2};[x,y]=retgranet(init,x*pas(1),y*pas(2));x=x/pas(1);y=y/pas(2);w{1}=x;w{2}=y;end;% transformation de granet 10 2011   
[a,d]=retcouche_cylindrique_radial(init,w,c);return;end;

if nargout==3;[a,d,w]=ret2couche(init,w,c,li,xsource);else [a,d]=ret2couche(init,w,c,li,xsource);end;
if imag(c)==-2; % stoquage dans la bibliotheque
a=retio(a,1);w=retio(w,1);
if isempty(bibli);bibli={iinit,{ww,real(c),a,w}};else bibli={bibli{:},{ww,real(c),a,w}};end;
end;
    

else                       %< ******************* 1D 
beta=init{2};si=init{3};cao=init{5};sog=init{end}.sog;n=size(beta,2);
if init{end}.granet==1;w{1}=retgranet(init,w{1});end;% transformation de granet   
x=w{1}/init{4};ep=w{2};
c=real(c);

pp=[];qq=[];q=[];
if length(x)<=1 & isempty(si) & isempty(cao); % <------- milieu homogene et pas de symetries et pas cao: diagonalisation analytique
d=retsqrt(-ep(1)/ep(3)+(beta.'.^2)./(ep(2)*ep(3)),1);
%d=retbidouille(d);% <<<<<<<<<<<<<<<<<< MODIF

p=eye(n);q=p/ep(3);pp=eye(n);qq=pp*ep(3);   
if c==1;
[fmu,fmmu,fep]=retc1(x,ep,beta);if size(w,2)==2;[a1,a2,k,aa2]=retc2(fmu,fmmu,fep,beta,cao);else;[a1,a2,k]=retc2(fmu,fmmu,fep,beta,cao,w{3},w{4});end;
a3=-inv(rettoeplitz(fmu))*diag(beta); % matrice a3 permettant le calcul de hy=(1/mu)*(du/dx)
if ~isempty(si);a3=a3*si{1};a1=a1*si{1};end
else
a3=[];a1=[];
end 

else                                     % <------- milieu non homogene ou symetries ou cao
  
[fmu,fmmu,fep]=retc1(x,ep,beta);if size(w,2)==2;[a1,a2,k,aa2]=retc2(fmu,fmmu,fep,beta,cao);else;[a1,a2,k]=retc2(fmu,fmmu,fep,beta,cao,w{3},w{4});end;
if ~isempty(si);          % < ............ utilisation des proprietes de symetrie
aa1=a1*si{1};a1=si{2}*aa1;a2=si{2}*a2*si{1};
end    
if length(x)<=1 & isempty(cao)  %<.....  milieu homogene dans le cas d'utilisation de symetries: diagonalisation analytique
d=retsqrt(diag(a2*a1),1);
%d=retbidouille(d);% <<<<<<<<<<<<<<<<<< MODIF

p=eye(init{1});q=a1;    
else;                           %<..... diagonalisation numerique    
[p,d]=retc3(a1,a2,c);
end;                            %<..... diagonalisation numerique ?
if c>=0;q=a1*p;pp=retpinv(p);qq=retpinv(q);end;

if c==1;
%a3=-inv(rettoeplitz(fmu))*diag(beta);% matrice a3 permettant le calcul de hy=i/mu*du/dx
a3=-inv(rettoeplitz(fmu))*k;% matrice a3 permettant le calcul de hy=i/mu*du/dx avec changement de coordonnees
if ~isempty(si);a3=a3*si{1};a1=aa1;end
else
a3=[];if ~isempty(si);a1=aa1;end
end;
end;                                    % <------- milieu homogene et pas de symetries   ? ?

a={p,pp,q,qq,d,a3,a1,w,init{end}.d,init{2},struct('dim',1,'genre',2,'type',1,'sog',sog)};
if nargout==3;w={w{1},w{2},a1,aa2};end;
end;                       %< ******************* 2D 1D 

end;           % < +++++++  metaux ou dielectriques

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [a,d,w]=ret2couche(init,w,c,li,xsource) 
io=imag(c)~=0;c=real(c);
beta=init{2};ordre=init{3};nx=init{6};ny=init{7};n1=nx*ny;si=init{8};pas=init{end}.d;n=init{1};cao=init{10};delt=init{end}.delt;sog=init{end}.sog;
x=w{1};y=w{2};u=w{3};
if init{end}.granet;[x,y]=retgranet(init,x(1:end-1)*pas(1),y(1:end-1)*pas(2));x=[x/pas(1),1];y=[y/pas(2),1];w{1}=x;w{2}=y;end;% transformation de granet
% attention aux erreurs d'arrondi sur x(end)=1   y(end)=1 
nw=(size(w,2)==3);  % on recalcule eepz... 


kbidouillemax=16;
kbidouille=1;
while (kbidouille>0)&(kbidouille<kbidouillemax);%  <-----------BOUCLE SUR KBIDOUILLE

if ~nw;             % on reutilise des valeurs deja calculees
eepz=w{4};mmuz=w{5};muy=w{6};mux=w{7};epy=w{8};epx=w{9};hx=w{10};ez=w{11};hz=w{12};fex=w{13};fhx=w{14};fey=w{15};fhy=w{16};
else;
    
    
eepz=[];mmuz=[];muy=[];mux=[];epy=[];epx=[];hx=[];ez=[];hz=[];fex=[];fhx=[];fey=[];fhy=[];
end;    

% parm=max([1,sqrt(max(sum(abs(beta.^2),1)))*.2]); % pour retbidouille
parm=kbidouille*max([1,sqrt(max(sum(abs(beta.^2),1)))*.2]); % pour retbidouille

%  < _________couche non homogene ou changement de coordonnees
if ~(all(all(all(repmat(u(1,1,:),size(u,1),size(u,2))==u)))&(u(1)==u(2))&(u(1)==u(3))&(u(4)==u(5))&(u(4)==u(6))&isempty(cao)); %guide

%calcul de la matrice m=[0,a1;a2,0] associee au passage dans un milieu periodique 2D invariant en z
if ~isempty(cao); % changement de coordonnees
kx=retmat(apodisecao(cao{1})*retdiag(beta(1,1:nx)),ny);
ky=retmat(apodisecao(cao{2})*retdiag(beta(2,1:nx:n1)),-nx);
else
kx=retdiag(beta(1,:));ky=retdiag(beta(2,:));
end;
[prv,version]=retversion;if version>7.003;nxnymax=inf;else;nxnymax=4000;end;   % limite du calcul complexe dans l'inversion de ret2li. Peut creer des instabilités
if ~isempty(si);  % < XXXX  utilisation des proprietes de symetrie 

if isempty(eepz);eepz=ret2li(x,y,1./u(:,:,6),nx,ny,0,0,li(1),si,nxnymax,3);else;eepz=retio(eepz);end;  %calcul de eepz
a1=i*(retprod(si{2}(:,1:n1)*kx,eepz,ky*si{3}(1:n1,:))-retprod(si{2}(:,1:n1)*kx,eepz,kx*si{3}(n1+1:2*n1,:))+retprod(si{2}(:,n1+1:2*n1)*ky,eepz,ky*si{3}(1:n1,:))-retprod(si{2}(:,n1+1:2*n1)*ky,eepz,kx*si{3}(n1+1:2*n1,:)));
if c==1 | li(7)==2 ;ez=retprod(eepz,ky*si{3}(1:nx*ny,:))-retprod(eepz,kx*si{3}(nx*ny+1:end,:));ez=retio(ez,io);end; %tableaux utilises pour calculer les champs
if (nargout==3)&nw;eepz=retio(eepz,io);else;retio(eepz,-1);clear eepz;end;

if isempty(muy);muy=ret2li(x,y,u(:,:,2),nx,ny,1,0,li(3),si,nxnymax,5);else;muy=retio(muy);end; %muy
a1=retfullsparse(a1);a1=a1+i*retprod(si{2}(:,1:n1),muy,si{3}(n1+1:2*n1,:));
if (nargout==3)&nw;muy=retio(muy,io);else;retio(muy,-1);clear muy;end;

if isempty(mux);mux=ret2li(x,y,u(:,:,1),nx,ny,0,1,li(4),si,nxnymax,4);else;mux=retio(mux);end;  % mux
a1=retfullsparse(a1);a1=a1-i*retprod(si{2}(:,n1+1:2*n1),mux,si{3}(1:n1,:));
a1=full(a1);
if (nargout==3)&nw;mux=retio(mux,io);else;retio(mux,-1);clear mux;end;

if isempty(mmuz);mmuz=ret2li(x,y,1./u(:,:,3),nx,ny,0,0,li(2),si,nxnymax,6);else;mmuz=retio(mmuz);end;  %calcul de mmuz
a2=i*(retprod(-si{4}(:,1:n1)*kx,mmuz,ky*si{1}(1:n1,:))+retprod(si{4}(:,1:n1)*kx,mmuz,kx*si{1}(n1+1:2*n1,:))-retprod(si{4}(:,n1+1:2*n1)*ky,mmuz,ky*si{1}(1:n1,:))+retprod(si{4}(:,n1+1:2*n1)*ky,mmuz,kx*si{1}(n1+1:2*n1,:)));

if c==1 | li(8)==2 ;hz=-retprod(mmuz,ky*si{1}(1:n1,:))+retprod(mmuz,kx*si{1}(n1+1:2*n1,:));hz=retio(hz,io);end; %tableaux utilises pour calculer les champs
if (nargout==3)&nw;mmuz=retio(mmuz,io);else;retio(mmuz,-1);clear mmuz;end;

if isempty(epy);epy=ret2li(x,y,u(:,:,5),nx,ny,1,0,li(5),si,nxnymax,2);else;epy=retio(epy);end; %calcul de epy
a2=retfullsparse(a2);a2=a2-i*retprod(si{4}(:,1:n1),epy,si{1}(n1+1:2*n1,:));
if (nargout==3)&nw;epy=retio(epy,io);else;retio(epy,-1);clear epy;end;
 
if isempty(epx);epx=ret2li(x,y,u(:,:,4),nx,ny,0,1,li(6),si,nxnymax,1);else;epx=retio(epx);end;  %calcul de epx
a2=retfullsparse(a2);a2=a2+i*retprod(si{4}(:,n1+1:2*n1),epx,si{1}(1:n1,:));
if (nargout==3)&nw;epx=retio(epx,io);else;retio(epx,-1);clear epx;end;clear kx ky;
a2=full(a2);


else;  % < XXXX  non utilisation des proprietes de symetrie 

if isempty(eepz);  %calcul de eepz
eepz=ret2li(x,y,1./u(:,:,6),nx,ny,0,0,li(1));else;eepz=retio(eepz);
end;
a1=1i*[[kx*eepz*ky,-kx*eepz*kx];[ky*eepz*ky,-ky*eepz*kx]];
if c==1 | li(7)==2 ;ez=[eepz*ky,-eepz*kx];ez=retio(ez,io);end; %tableaux utilises pour calculer les champs

if (nargout==3)&nw;eepz=retio(eepz,io);else;clear eepz;end;

if isempty(muy);muy=ret2li(x,y,u(:,:,2),nx,ny,1,0,li(3));else;muy=retio(muy);end; %muy
a1=retfullsparse(a1);a1(1:n1,n1+1:2*n1)=a1(1:n1,n1+1:2*n1)+i*muy;
if (nargout==3)&nw;muy=retio(muy,io);else;clear muy;end;

if isempty(mux);mux=ret2li(x,y,u(:,:,1),nx,ny,0,1,li(4));else;mux=retio(mux);end;  % mux
a1=retfullsparse(a1);a1(n1+1:2*n1,1:n1)=a1(n1+1:2*n1,1:n1)-i*mux;
if (nargout==3)&nw;mux=retio(mux,io);else;clear mux;end;

if isempty(mmuz);mmuz=ret2li(x,y,1./u(:,:,3),nx,ny,0,0,li(2));else;mmuz=retio(mmuz);end;  %calcul de mmuz
a2=1i*[[-kx*mmuz*ky,kx*mmuz*kx];[-ky*mmuz*ky,ky*mmuz*kx]];
if c==1 | li(8)==2 ;hz=[-mmuz*ky,mmuz*kx];hz=retio(hz,io);end; %tableaux utilises pour calculer les champs
if (nargout==3)&nw;mmuz=retio(mmuz,io);else;clear mmuz;end;

if isempty(epy);epy=ret2li(x,y,u(:,:,5),nx,ny,1,0,li(5));else;epy=retio(epy);end; %calcul de epy
a2=retfullsparse(a2);a2(1:n1,n1+1:2*n1)=a2(1:n1,n1+1:2*n1)-i*epy;
if (nargout==3)&nw;epy=retio(epy,io);else;clear epy;end;
 
if isempty(epx);epx=ret2li(x,y,u(:,:,4),nx,ny,0,1,li(6));else;epx=retio(epx);end;  %calcul de epx
a2=retfullsparse(a2);a2(n1+1:2*n1,1:n1)=a2(n1+1:2*n1,1:n1)+i*epx;
if (nargout==3)&nw;epx=retio(epx,io);else;clear epx,end;clear kx ky;
    
    
end;   % < XXXX   utilisation des proprietes de symetrie  ?
%a1,a2
moharam=0;sens=0;
if isempty(cao);% singularites ???????????????????????????????????????????????????????????????????????
if ~isempty(si);ne=size(si{6},1);nh=size(si{8},1);else;ne=init{1}/2;nh=ne;end
if size(u,1)==1 & size(u,2)>1  & nx==1 & ne>0;moharam=nh;sens=-1;end;% 
if size(u,2)==1 & size(u,1)>1   & ny==1 & nh>0;moharam=ne;sens=1;end;
end;
if c~=-2;[p,d]=retc3(a1,a2,c,moharam,sens);if c~=-3;q=a1*p;else;q=[];end;else p=[];q=[];d=[];end;  % diagonalisation

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   utile pour les milieux homogenes non isotropes ...
%        % orthogonalisation au sens du flux de Poynting des modes de vp egales
%        bb=find(abs(real(d))<1.e-6&abs(imag(d))>eps);
%        if ~isempty(bb);
%        q=a1*p;
%        if isempty(si);p1=p(:,bb);q1=q(:,bb);else;p1=si{3}*p(:,bb);q1=si{1}*q(:,bb);end;
%        n2=length(bb);masque=zeros(1,n2);
%        for ii=1:n2;
%        if masque(ii)==0;
%        jj=find(d(bb)==d(bb(ii)));masque(jj)=1;    
%        if length(jj)>1;
%        b=bb(jj);    
%        [cho,test]=chol(p1(:,jj)'*[[zeros(n1/2),-eye(n1/2)];[eye(n1/2),zeros(n1/2)]]*q1(:,jj)*diag(d(b)));
%        if test==0;p(:,b)=p(:,b)*inv(cho);end;
%        end;end;end;
%        end;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


else  % < _________  espace homogene et isotrope

if c==1; %< ....tableaux utilises pour calculer les champs Ez et Hz
kx=retdiag(beta(1,:));ky=retdiag(beta(2,:));
if isempty(eepz);eepz=ret2li(x,y,1./u(:,:,6),nx,ny,0,0,li(7));else;eepz=retio(eepz);end;
if isempty(mmuz);mmuz=ret2li(x,y,1./u(:,:,3),nx,ny,0,0,li(8));else;mmuz=retio(mmuz);end;
ez=[eepz*ky,-eepz*kx];if ~isempty(si);ez=ez*si{3};end;ez=retio(ez,io);
hz=[-mmuz*ky,mmuz*kx];if ~isempty(si);hz=hz*si{1};end;hz=retio(hz,io);
end;     %< ....
if (nargout==3) & nw;% <---
muy=ret2li(x,y,u(:,:,2),nx,ny,1,0,li(3));muy=retio(muy,io);  % muy
mux=ret2li(x,y,u(:,:,1),nx,ny,0,1,li(4));mux=retio(mux,io);  % mux
epy=ret2li(x,y,u(:,:,5),nx,ny,1,0,li(5));epy=retio(epy,io);  % epy
epx=ret2li(x,y,u(:,:,4),nx,ny,0,1,li(6));epx=retio(epx,io);  % epx
eepz=ret2li(x,y,1./u(:,:,6),nx,ny,0,0,li(7));eepz=retio(eepz,io);
mmuz=ret2li(x,y,1./u(:,:,3),nx,ny,0,0,li(8));mmuz=retio(mmuz,io);
else;              % <---
retio(eepz,-1);retio(mux,-1);retio(muy,-1);retio(mmuz,-1);retio(epx,-1);retio(epy,-1);    
clear eepz mux muy mmuz epx epy
end;               % <---
clear kx ky;



if isempty(si); % < XXXX  pas de symetrie   decomposition E et H pas de symetrie
ux=beta(1,:);uy=beta(2,:);
d=repmat(retsqrt(ux.^2+uy.^2-u(1,1,1)*u(1,1,4),1).',2,1);
%d=retbidouille(d);% <<<<<<<<<<<<<<<<<< MODIF

uu=sqrt(ux.^2+uy.^2);
f=find(abs(uu)<100*eps);ux(f)=cos(delt);uy(f)=sin(delt);uu(f)=1;% onde normale  cas degenere angle delt par defaut   
ux=ux./uu;uy=uy./uu;
% p=[[diag(ux);diag(uy)],[diag(uy);diag(-ux)]];% modif 12 2011
d2=-retbidouille(d([1:n/2]).',1,parm).^2;

p=[[diag(ux);diag(uy)],[diag(uy);diag(-ux)]];
q=i*[[diag(uy);diag(-ux)]*u(1,1,1),[diag(-ux);diag(-uy)]*diag(d2/u(1,1,4))];

else;  % < XXXX   cas des symetries
p=eye(n);q=zeros(n);p1=si{3};d=zeros(n,1);
masque=zeros(1,n);
for ii=1:n;  % <*** boucle ii
if masque(ii)==0;  % <+++    masque=0
masque(ii)=1;bm=find(masque==0);
ff=find(p1(:,ii)~=0);
f=ff(1);if f>n1;f=f-n1;end;
k=find((p1(f,bm)~=0)|(p1(f+n1,bm)~=0));
masque(bm(k))=1;
ux=beta(1,f);uy=beta(2,f); 
d2=retsqrt(ux^2+uy^2-u(1,1,1)*u(1,1,4),1);
d([ii,bm(k)])=d2;
%d([ii,bm(k)])=retbidouille(d2);% <<<<<<<<<<<<<<<<<< MODIF
d2=-retbidouille(d2,1,parm).^2;

uu=sqrt(ux^2+uy^2);
if abs(uu)<100*eps; ux=cos(delt);uy=sin(delt);else ux=ux/uu;uy=uy/uu;end; % onde normale  cas degenere angle delt par defaut   
 
if ~isempty(k);      % <---- k~=[]
kk=bm(k(1));    
prv=(p1([f,f+n1],[ii,kk])\([[ux;uy],[uy;-ux]]));  %passage E//  H//
p(:,[ii,kk])=p(:,[ii,kk])*prv;
q(:,[ii,kk])=si{2}*[p1(1+n1:2*n1,[ii,kk]);-p1(1:n1,[ii,kk])]*prv*[[i*u(1,1,1);0],[0;i*d2/u(1,1,4)]];
else          % <---une seule op dans le mode
if abs(ux*p1(f,ii)+uy*p1(f+n1,ii))<abs(uy*p1(f,ii)-ux*p1(f+n1,ii));%cas H//
prv=i*d2/u(1,1,4);
else;%cas E//
prv=i*u(1,1,1);
end;
q(:,ii)=prv*si{2}*([p1(1+n1:2*n1,ii);-p1(1:n1,ii)]);
end;                  % <---- k~=[] ?
end;  % <+++    masque=0

end;    % <*** boucle ii
p=sparse(p);q=sparse(q);
end;% < XXXX  symetries ?
end; % < _________

rc=min(rcond(full(p)),rcond(full(q)));if rc<=5*eps;kbidouille=2*kbidouille;else kbidouille=0;end;

end;  %  <-----------BOUCLE SUR KBIDOUILLE (on recommence si p ou q mal conditionnees)



if c<0;pp=[];qq=[];else pp=retpinv(p);qq=retpinv(q);end;

if c==1; %tableaux utilises pour calculer les champs non encore calcules
mx=size(x,2);my=size(y,2);mxx=1;myy=1;
if ~isempty(xsource);  % tableaux reduits pour les sources 
xxx=mod(xsource(1)/pas(1),1);yyy=mod(xsource(2)/pas(2),1);

for iii=1:mx;if iii==1 x0=0;else x0=x(iii-1);end;
fx=find((xxx>=x0)&(xxx<=x(iii)));    
for jjj=1:my;if jjj==1 y0=0;else y0=y(jjj-1);end;
fy=find((yyy>=y0)&(yyy<=y(jjj)));
if ~isempty(fx)&~isempty(fy) iiii=iii;jjjj=jjj;end;
end;end;
else;mxx=1:mx;myy=1:my;iiii=mxx;jjjj=myy;% tableaux complets
end;  %tableaux utilises pour calculer les champs non encore calcules

if isempty(fex);
l_mxx=length(mxx);l_myy=length(myy);
fex=cell(1,l_myy);
for ii=1:l_myy; % fonctions discontinues en x devenant continues
trouve=1;for jj=1:ii-1;if all(u(:,jjjj(myy(ii)),4)==u(:,jjjj(myy(jj)),4));fex{myy(ii)}=fex{myy(jj)};trouve=0;break;end;end;	
if trouve;fex{myy(ii)}=retio(ret2li(x,1,u(:,jjjj(myy(ii)),4),nx,ny,0,1,li(9)),io);end;
end;
fey=cell(1,l_mxx);
for ii=1:l_mxx; % fonctions discontinues en y devenant continues 
trouve=1;for jj=1:ii-1;if all(u(iiii(mxx(ii)),:,5)==u(iiii(mxx(jj)),:,5));fey{mxx(ii)}=fey{mxx(jj)};trouve=0;break;end;end;	
if trouve;fey{mxx(ii)}=retio(ret2li(1,y,u(iiii(mxx(ii)),:,5),nx,ny,1,0,li(11)),io);end;
end;
fhx=cell(1,l_myy);
for ii=1:l_myy; % fonctions discontinues en x devenant continues
trouve=1;for jj=1:ii-1;if all(u(:,jjjj(myy(ii)),1)==u(:,jjjj(myy(jj)),1));fhx{myy(ii)}=fhx{myy(jj)};trouve=0;break;end;end;	
if trouve;fhx{myy(ii)}=retio(ret2li(x,1,u(:,jjjj(myy(ii)),1),nx,ny,0,1,li(10)),io);end;
end;
fhy=cell(1,l_mxx);
for ii=1:l_mxx; % fonctions discontinues en y devenant continues
trouve=1;for jj=1:ii-1;if all(u(iiii(mxx(ii)),:,2)==u(iiii(mxx(jj)),:,2));fhy{mxx(ii)}=fhy{mxx(jj)};trouve=0;break;end;end;	
if trouve;fhy{mxx(ii)}=retio(ret2li(1,y,u(iiii(mxx(ii)),:,2),nx,ny,1,0,li(12)),io);end;
end;
end;
end;

a={p,pp,q,qq,d,ez,hz,fex,fhx,fey,fhy,w,pas,init{2},init{11},struct('dim',2,'genre',2,'type',1,'sog',sog)};


%end;  %  <-----------BOUCLE SUR KBIDOUILLE (on recommence si p ou q mal conditionnees)% (pourquoi j'avais mis cette ligne la ???????????????)

if (nargout==3)& nw;
w={x/pas(1),y/pas(2),w{3},eepz,mmuz,muy,mux,epy,epx,hx,ez,hz,fex,fhx,fey,fhy};% attention si transformee de granet on remet dans w ev sortie les x y physiques !
                                                                       % Le w dans a a les coordonnees numeriques il peut etre utilise par retchamp
end;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [aa,dd]=retmcouche(init,www,cc);
%  function [aa,dd]=retmcouche(init,www,cc);
%  calcul du descripteur de texture aa utilisee par  retc pour le metal infiniment conducteur

n=init{1};beta=init{2};x=www{1};dx=www{2};ep=www{3};pol=www{4};mm=www{5};

if isempty(x);dd=[];% metal non troue
if pol==0; % metal electrique
%g1=[eye(n),zeros(n)];g2=zeros(n,0);
g1=[eye(n),1.e-5*eps*eye(n)];g2=zeros(n,0);% modif 10 2010 pour la stabilite
else; % metal magnetique
%g1=[zeros(n),eye(n)];g2=zeros(n,0);
g1=[1.e-5*eps*eye(n),eye(n)];g2=zeros(n,0);% modif 10 2010 pour la stabilite
end;    
aa={g1,g2,[],0,n,cell(1,5),www,struct('dim',init{end}.dim,'genre',2,'type',2,'sog',0)};
return;
end 


nt=size(x,1);% nombre de trous
champ=cell(1,5);
if init{end}.dim==1; %1D
               %%%%%%%%%%%%%%%%%%%%
               %      1  D        % 
               %%%%%%%%%%%%%%%%%%%%

si=init{3};d=init{end}.d;
n1=size(beta,2);if isempty(mm);mm=repmat(n1,nt,1);end;
mmax=sum((mm+pol/2));
sc=zeros(mmax,n1);scm=zeros(mmax,n1);
dd=zeros(mmax,1);ww=zeros(mmax,1);
if cc==1;cha=zeros(mmax,4,2);chx=zeros(mmax,1);chdx=zeros(mmax,1);chalp=zeros(mmax,1);chtp=zeros(mmax,1);end;

% calcul de sc et sm

parm=pol;tp=1;
m=0;w=[];  %  nombre de modes
for it=1:nt;
for mx=1-pol/2:mm(it);
m=m+1;
alp=mx*pi/dx(it);bet=beta;x0=x(it);d0=dx(it);
[s,c,sm,cm]=calsc(d0,bet,x0,dx,alp,tp,parm);
sc(m,:)=s;scm(m,:)=sm;
dd(m)=retbidouille(retsqrt(alp^2-ep(it,1)/ep(it,3),1),1);
w(m)=dd(m)*ep(it,3);

if cc==1;
% amplitudes des champs: 1 montants: 2 descendants
cha(m,:,1)=[1,i*w(m),i*alp*ep(it,3),i*dd(m)];
cha(m,:,2)=[1,-i*w(m),i*alp*ep(it,3),-i*dd(m)];
chx(m)=x(it);     % centre
chdx(m)=dx(it);    % largeur
chalp(m)=alp;     % alpha
chtp(m)=pol-1;% tp: 1 cos  -1 sin
end;

end;
end;
sc=sc(1:m,:);scm=scm(1:m,:);dd=dd(1:m);w=w(1:m);
if cc==1;cha=cha(1:m,:,:);chx=chx(1:m);chdx=chdx(1:m);chalp=chalp(1:m);chtp=chtp(1:m);end;
% construction de la matrice G

if pol==0;
g1=diag(1./w)*sc;g2=scm.'/d;
if ~isempty(si);g1=g1*si{1};g2=si{2}*g2;end; %  utilisation des symetries
f=find(max(abs(g2).^2,[],1)>eps);m=size(f,2);% elimination des modes nulls (dans le cas des symetries)
g1=g1(f,:);g2=g2(:,f);dd=dd(f);
g1=[[zeros(m,n),g1];[eye(n),zeros(n)]];g2=[[-eye(m),eye(m)];[g2,g2]];
else; 
g1=sc;g2=scm.'*diag(w)/d;
if ~isempty(si);g1=g1*si{1};g2=si{2}*g2;end; %  utilisation des symetries
f=find(max(abs(g1).^2,[],2)>eps);m=size(f,1);% elimination des modes nulls (dans le cas des symetries)
g1=g1(f,:);g2=g2(:,f);dd=dd(f);
g1=[[g1,zeros(m,n)];[zeros(n),eye(n)]];g2=[[eye(m),eye(m)];[-g2,g2]];
end;
if cc==1;cha=cha(f,:,:);chx=chx(f);chdx=chdx(f);chalp=chalp(f);chtp=chtp(f);end;
else;
               %%%%%%%%%%%%%%%%%%%%
               %      2  D        % 
               %%%%%%%%%%%%%%%%%%%%
nx=init{6};ny=init{7};si=init{8};d=init{end}.d;
n1=nx*ny;if isempty(mm);mm=repmat([nx,ny],nt,1);end;
sc=zeros(0,n1);cs=zeros(0,n1);scm=zeros(0,n1);csm=zeros(0,n1);
alpha=zeros(0,2);
u=zeros(0,1);v=zeros(0,1);w=zeros(0,1);dd=zeros(0,1);
if cc==1;chax1=zeros(0,6);chay1=zeros(0,6);chax2=zeros(0,6);chay2=zeros(0,6);chx=zeros(0,2);chdx=zeros(0,2);chalp=zeros(0,2);chtp=zeros(0,2);end;

% calcul de sc scm cs et csm (beta peut etre complexe donc scm n'est pas toujours conj(sc) ...)
parm=1;
m=0;  %  nombre de modes
for it=1:nt;
if dx(it,1)>d(1);mmx=nx-1;parmx=0;else;mmx=mm(it,1);parmx=1;end; 
if dx(it,2)>d(2);mmy=ny-1;parmy=0;else;mmy=mm(it,2);parmy=1;end; 
for mx=0:mmx;

if parmx==1;
alp=mx*pi/dx(it,1);tp=1;d0=dx(it,1);
else;
alp=beta(1,mx+1);tp=2;d0=d(1);    
end;
bet=beta(1,1:nx);x0=x(it,1);[s,c,sm,cm]=calsc(d0,bet,x0,dx,alp,tp,parm);

sx=s;cx=c;sxm=sm;cxm=cm;alpx=alp;

for my=0:mmy;
m=m+1; % creation d'un mode
if parmy==1;
alp= my*pi/dx(it,2);tp=1;d0=dx(it,2);
else;
alp=beta(2,1+nx*my);tp=2;d0=d(2);
end;
bet=beta(2,1:nx:n1);x0=x(it,2);[s,c,sm,cm]=calsc(d0,bet,x0,dx,alp,tp,parm);
alpha=[alpha;[alpx,alp]];

sc=[sc;reshape(sx.'*c,1,n1)];cs=[cs;reshape(cx.'*s,1,n1)];
scm=[scm;reshape(sxm.'*cm,1,n1)];csm=[csm;reshape(cxm.'*sm,1,n1)];

n2=ep(it,3)*ep(it,6);ww=i/ep(it,6-3*pol/2);
dd=[dd;retbidouille(retsqrt(alpha(m,1)^2+alpha(m,2)^2-n2,1),1)];
u=[u;ww*alpha(m,1)*alpha(m,2)/dd(m)];
v=[v;ww*(alpha(m,2)^2-n2)/dd(m)];
w=[w;ww*(alpha(m,1)^2-n2)/dd(m)];

if cc==1;
% amplitudes des champs
if pol==0; % metal electrique
chax1=[chax1;[u(m),v(m),-i/ep(it,6)*alpha(m,2),-1,0,-alpha(m,1)/dd(m)]];
chay1=[chay1;[-w(m),-u(m),i/ep(it,6)*alpha(m,1),0,-1,-alpha(m,2)/dd(m)]];
chax2=[chax2;[u(m),v(m),i/ep(it,6)*alpha(m,2),1,0,-alpha(m,1)/dd(m)]];
chay2=[chay2;[-w(m),-u(m),-i/ep(it,6)*alpha(m,1),0,1,-alpha(m,2)/dd(m)]];
else; % metal magnetique
chax1=[chax1;[1,0,alpha(m,1)/dd(m),u(m),v(m),-i/ep(it,3)*alpha(m,2)]];
chay1=[chay1;[0,1,alpha(m,2)/dd(m),-w(m),-u(m),i/ep(it,3)*alpha(m,1)]];
chax2=[chax2;[1,0,-alpha(m,1)/dd(m),-u(m),-v(m),-i/ep(it,3)*alpha(m,2)]];
chay2=[chay2;[0,1,-alpha(m,2)/dd(m),w(m),u(m),i/ep(it,3)*alpha(m,1)]];
end;
chx=[chx;x(it,:)];          % centre
chdx=[chdx;min(dx(it,:),d)];           % largeur
chalp=[chalp;alpha(m,:)];     % alpha
chtp=[chtp;[(2-parmx)*(1-pol),(2-parmy)*(1-pol)]];%  si pol=0 (electrique): tp= 2 (exp)  tp=  1(cos) 
                                                  %  si pol=2 (magnetique): tp=-2 (exp)  tp= -1(cos)
end;

end;end;end;
if cc==1;cha=cat(3,[chax1;chay1],[chax2;chay2]);chx=[chx;chx];chdx=[chdx;chdx];chalp=[chalp;chalp];chtp=[chtp;chtp];end

% construction de la matrice G
g1=[[sc;zeros(m,n1)],[zeros(m,n1);cs]];
g2=[[csm.'*diag(u);scm.'*diag(v)],[-csm.'*diag(w);-scm.'*diag(u)]];
if pol==0; % metal electrique
if ~isempty(si);g1=g1*si{3};g2=si{2}*g2;end; %  utilisation des symetries
f=find(max(abs(g1).^2,[],2)>eps);m=size(f,1);% elimination des modes nulls 
dd=[dd;dd];dd=dd(f);
g1=g1(f,:);g2=g2(:,f)/(d(1)*d(2));
g1=[[zeros(m,n),g1];[eye(n),zeros(n)]];g2=[[-eye(m),eye(m)];[g2,g2]];
else; % metal magnetique
if ~isempty(si);g1=g1*si{1};g2=si{4}*g2;end; %  utilisation des symetries
f=find(max(abs(g1).^2,[],2)>eps);m=size(f,1);% elimination des modes nulls 
dd=[dd;dd];dd=dd(f);
g1=g1(f,:);g2=g2(:,f)/(d(1)*d(2));
g1=[[g1,zeros(m,n)];[zeros(n),eye(n)]];g2=[[eye(m),eye(m)];[g2,-g2]];
end;
if cc==1;cha=cha(f,:,:);chx=chx(f,:);chdx=chdx(f,:);chalp=chalp(f,:);chtp=chtp(f,:);end;
end;
if cc==1;champ={chx,chdx,chalp,chtp,cha};end;
aa={sparse(g1),sparse(g2),dd,m,n,champ,www,struct('dim',init{end}.dim,'genre',2,'type',2,'sog',0)};



function [s,c,sm,cm]=calsc(d0,bet,x0,dx,alp,tp,parm);% calcul de 'l'integrale sin*epp..' si parm=0 s sm  (E//)  si parm=2 c cm (H//)  s1 parm=1  les 2 (2D)
c=[];cm=[];
if tp==2;c=(bet==alp)*sqrt(d0);s=i*c;cm=c;sm=-s;else;
prv1=retsinc((bet+alp)*d0/(2*pi));prv2=retsinc((bet-alp)*d0/(2*pi));
if parm==1;
c=sqrt(d0/2)*exp(i*bet*x0).*(exp(i*alp*d0/2)*prv1+exp(-i*alp*d0/2)*prv2);
cm=sqrt(d0/2)*exp(-i*bet*x0).*(exp(i*alp*d0/2)*prv2+exp(-i*alp*d0/2)*prv1);
s=-i*sqrt(d0/2)*exp(i*bet*x0).*(exp(i*alp*d0/2)*prv1-exp(-i*alp*d0/2)*prv2);
sm=-i*sqrt(d0/2)*exp(-i*bet*x0).*(exp(i*alp*d0/2)*prv2-exp(-i*alp*d0/2)*prv1);
if alp==0;c=c/sqrt(2);cm=cm/sqrt(2);end;
end;
if parm==2;
s=sqrt(d0/2)*exp(i*bet*x0).*(exp(i*alp*d0/2)*prv1+exp(-i*alp*d0/2)*prv2);
sm=sqrt(d0/2)*exp(-i*bet*x0).*(exp(i*alp*d0/2)*prv2+exp(-i*alp*d0/2)*prv1);
if alp==0;s=s/sqrt(2);sm=sm/sqrt(2);end;
end;
if parm==0;
s=-i*sqrt(d0/2)*exp(i*bet*x0).*(exp(i*alp*d0/2)*prv1-exp(-i*alp*d0/2)*prv2);
sm=-i*sqrt(d0/2)*exp(-i*bet*x0).*(exp(i*alp*d0/2)*prv2-exp(-i*alp*d0/2)*prv1);
end;
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function ffff=ret2li(xx,yy,u,nx,ny,cx,cy,sens,si,nxnymax,pol);
% calcul de la matrice ffff
%  tf(g)=ffff*tf(f)
%     g=u*f
% cx:f continue en x  cy:f continue en y
% les continuitees de g sont inversees par multiplication par u
% sens:parametre indiquant le sens du calcul(facultatif)
%

if nargin<8 sens=1;end;
if nargin<9 si=[];end;
if nargin<10 nxnymax=inf;end;
if sens==2;cx=1;cy=1;sens=1;end;
% essai de simplification de xx yy u
mx=size(xx,2);my=size(yy,2);
x=xx(mx);uu=u(mx,:);
for ii=mx-1:-1:1;if ~all(u(ii,:)==u(ii+1,:));uu=[u(ii,:);uu];x=[xx(ii),x];end;end;
u=uu(:,my);y=yy(my);
for jj=my-1:-1:1;if ~all(uu(:,jj)==uu(:,jj+1));u=[uu(:,jj),u];y=[yy(jj),y];end;end;

n=nx*ny;
if length(retelimine(u(:)))==1;ffff=u(1)*speye(n);return;end;% cas constant (mu sans pml)
mx=size(x,2);my=size(y,2);	
alx=(-nx+1:nx-1)*(2*pi);aly=(-ny+1:ny-1)*(2*pi);
if nargin<8;sens=0;end;
if sens==0;sens=cx|(~cy);end;% permet d'eviter l' inversion finale (sauf si cx=cy=0)

topx=rettoeplitz(1:2*nx-1);topy=rettoeplitz(1:2*ny-1);% pour gain de temps

if sens==1 % de X en Y
if cy;f=retf(y,u,aly);else f=retf(y,1./u,aly);end;    
ff=zeros(ny,ny,mx);
if (cx&cy)|((~cx)&(~cy)) 
for ii=1:mx;ff(:,:,ii)=rettoeplitz(f(ii,:),topy);end;%test1=retcompare(ff,rettoeplitz(f,topy,[]))
else
for ii=1:mx;ff(:,:,ii)=inv(rettoeplitz(f(ii,:),topy));end;
end;
fff=reshape(retf(x,reshape(ff,ny^2,mx),alx),ny,ny,2*nx-1);

if n<nxnymax;  % <-------version 'normale' 
% ffff=zeros(nx,nx,ny,ny);
% for ii=1:ny;for jj=1:ny; 
% ffff(:,:,ii,jj)=rettoeplitz(fff(ii,jj,:),topx);
% end;end;
% ffff=permute(ffff,[1,3,2,4]);
ffff=rettoeplitz(fff,topx,[1,3,2,4]);% modif 3 2014

ffff=reshape(ffff,n,n);
if ~cx;if isempty(si);ffff=inv(ffff);else;ffff=retblkinv(ffff,si{11},si{9},si{10},si{12}(pol));end;end;
% if ~cx; ffff=inv(ffff);end;

ffffi=[];

else;         % <----version 'economique' en memoire
% ffff=zeros(nx,nx,ny,ny);
% for ii=1:ny;for jj=1:ny;
% ffff(:,:,ii,jj)=rettoeplitz(real(fff(ii,jj,:)),topx);
% end;end
% ffff=permute(ffff,[1,3,2,4]);
ffff=rettoeplitz(real(fff),topx,[1,3,2,4]);  % modif 3 2014

ffff=reshape(ffff,n,n);ffff=retio(ffff,1);
% ffffi=zeros(nx,nx,ny,ny);
% for ii=1:ny;for jj=1:ny;
% ffffi(:,:,ii,jj)=rettoeplitz(imag(fff(ii,jj,:)),topx);
% end;end
% ffffi=permute(ffffi,[1,3,2,4]);
ffffi=rettoeplitz(imag(fff),topx,[1,3,2,4]);  % modif 3 2014

ffffi=reshape(ffffi,n,n);ffffi=retio(ffffi,1);
if ~cx;if isempty(si);[ffff,ffffi]=retinv(ffff,ffffi);else;[ffff,ffffi]=retblkinv(ffff,si{11},si{9},si{10},si{12}(pol),ffffi);end;end;
end;              % <-----version 'economique'  ?   
       
    
else  % de Y en X
ff=zeros(nx,nx,my);
if cx;f=retf(x,u.',alx);else;f=retf(x,1./u.',alx);end;
if (cx&cy)|((~cx)&(~cy));
for ii=1:my;ff(:,:,ii)=rettoeplitz(f(ii,:),topx);end;%test3=retcompare(ff,rettoeplitz(f,topx,[]))
else
for ii=1:my;ff(:,:,ii)=inv(rettoeplitz(f(ii,:),topx));end;
end;
fff=reshape(retf(y,reshape(ff,nx^2,my),aly),nx,nx,2*ny-1);

if n<nxnymax; % <----- version 'normale'
% ffff=zeros(ny,ny,nx,nx);
% for ii=1:nx;for jj=1:nx;
% ffff(:,:,ii,jj)=rettoeplitz(fff(ii,jj,:),topy);
% end;end;
% ffff=permute(ffff,[3,1,4,2]);
ffff=rettoeplitz(fff,topy,[3,1,4,2]);% modif 3 2014

ffff=reshape(ffff,n,n);
if ~cy;if isempty(si);ffff=inv(ffff);else;ffff=retblkinv(ffff,si{11},si{9},si{10},si{12}(pol));end;end; 
ffffi=[];

else;          % <----- version 'economique' en memoire
% ffff=zeros(ny,ny,nx,nx);
% for ii=1:nx;for jj=1:nx;
% ffff(:,:,ii,jj)=rettoeplitz(real(fff(ii,jj,:)),topy);
% %ffff(ii,:,jj,:)=rettoeplitz(real(fff(ii,jj,:)),topy);
% end;end;
% ffff=permute(ffff,[3,1,4,2]);
ffff=rettoeplitz(real(fff),topy,[3,1,4,2]);% modif 3 2014

ffff=reshape(ffff,n,n);ffff=retio(ffff,1);
% ffffi=zeros(ny,ny,nx,nx);
% for ii=1:nx;for jj=1:nx;
% ffffi(:,:,ii,jj)=rettoeplitz(imag(fff(ii,jj,:)),topy);
% %ffffi(ii,:,jj,:)=rettoeplitz(imag(fff(ii,jj,:)),topy);
% end;end;
% ffffi=permute(ffffi,[3,1,4,2]);
ffffi=rettoeplitz(imag(fff),topy,[3,1,4,2]);% modif 3 2014

ffffi=reshape(ffffi,n,n);ffffi=retio(ffffi,1);
%if ~cy;[ffff,ffffi]=retinv(ffff,ffffi);end;
if ~cy;if isempty(si);[ffff,ffffi]=retinv(ffff,ffffi);else;[ffff,ffffi]=retblkinv(ffff,si{11},si{9},si{10},si{12}(pol),ffffi);end;end;
end;
end;          % <-----version 'economique'  ?  


%if n>=nxnymax;ffff={ffff,ffffi};end;% ???????????
if n>=nxnymax;ffff=retio(ffff)+i*retio(ffffi);end;
ffff=retsparse(ffff);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function ffff=ret2li_matrice(xx,yy,u,nx,ny,mx,my,sens);
% u=ret2li_matrice(x,y,U,nx,ny,mx,my,sens)
	% calcul de la matrice ffff
	%  ffff*{tf(g)}={tf(f)}
	% si sens=1: {f} tableau de p fonctions continues en x qui multiplies par une matrice carree u devient
	% un tableau de fonctions continues en y
	% si sens=- 1: {f} tableau de p fonctions continues en y qui multiplies par une matrice carree u devient
	% un tableau de fonctions continues en x
% La matrice est decrite par x y {U} points de discontinuite en x et y et valeurs à droite
% U est un cell array: size(U)=[p,p]  (ici p=3 )
% on obtient la matrice INVERSE
%  G  = U * F  espace reel
%  u *  g  =  f  espace de fourier (matrices)
% sens = -1  F continu en x  G continue en y
% sens = 1 F continu en y  G continue en x

if nargin<8;sens=-1;end;
if iscell(xx);x=xx{1};else;x=xx;end;if iscell(yy);y=yy{1};else;y=yy;end;% pour intégration de gauss
alx=(-nx+1:nx-1)*(2*pi);aly=(-ny+1:ny-1)*(2*pi);
p=size(u,2);
u_elimine=cell(size(u));n_elimine=zeros(size(u));for iii=1:length(u_elimine(:));u_elimine{iii}=retelimine(u{iii});n_elimine(iii)=length(u_elimine{iii})==1;end;	

if all(n_elimine(:));% cte pour gagner du temps;
for iii=1:numel(n_elimine);n_elimine(iii)=u_elimine{iii}(1);end; 
if sens~=0;n_elimine=inv(n_elimine);end;	
ffff=cell(size(u));[ffff{:}]=deal(sparse(nx*ny,nx*ny));for iii=1:numel(n_elimine);ffff{iii}=n_elimine(iii)*speye(nx*ny);end;	
return;end;

topx=rettoeplitz(1:2*nx-1);topy=rettoeplitz(1:2*ny-1);% pour gain de temps
switch sens;
case 1;% <*******   ffff * Tf(fonction continue en y)= Tf(fonction continue en x)
ff=zeros(p*nx,p*nx,my);
for iii=1:p;for jjj=1:p;     % tf en x fonction de y
if 	n_elimine(iii,jjj);for iy=1:my;ff([1:nx]+(iii-1)*nx,[1:nx]+(jjj-1)*nx,iy)=u_elimine{iii,jjj}(1)*speye(nx);end;else;
f=retf(xx,u{iii,jjj}.',alx);
for iy=1:my;ff([1:nx]+(iii-1)*nx,[1:nx]+(jjj-1)*nx,iy)=rettoeplitz(f(iy,:),topx);end;
end;end;end;                      % iii jjj
for iy=1:my;ff(:,:,iy)=inv(ff(:,:,iy));end;
fff=cell(p);[fff{:}]=deal(zeros(nx,nx,ny));
for iii=1:p;for jjj=1:p;fff{iii,jjj}=ff([1:nx]+(iii-1)*nx,[1:nx]+(jjj-1)*nx,:);end;end;
for iii=1:p;for jjj=1:p;fff{iii,jjj}=reshape(retf(yy,reshape(fff{iii,jjj},nx^2,my),aly),nx,nx,2*ny-1);end;end;
ffff=cell(p);[ffff{:}]=deal(zeros(ny,ny,nx,nx));
for iii=1:p;for jjj=1:p;
% for ii=1:nx;for jj=1:nx;ffff{iii,jjj}(:,:,ii,jj)=rettoeplitz(fff{iii,jjj}(ii,jj,:),topy);end;end;
% ffff{iii,jjj}=permute(ffff{iii,jjj},[3,1,4,2]);	
ffff{iii,jjj}=rettoeplitz(fff{iii,jjj},topy,[3,1,4,2]);% modif 3 2014
%for ii=1:nx;for jj=1:nx;ffff{iii,jjj}(ii,:,jj,:)=rettoeplitz(fff{iii,jjj}(ii,jj,:),topy);end;end;

ffff{iii,jjj}=reshape(ffff{iii,jjj},nx*ny,nx*ny);
end;end;

case -1     % <*******   ffff * Tf(fonction continue en x)= Tf(fonction continue en y)
ff=zeros(p*ny,p*ny,mx);
for iii=1:p;for jjj=1:p; % tf en y fonction de x
if 	n_elimine(iii,jjj);for ix=1:mx;ff([1:ny]+(iii-1)*ny,[1:ny]+(jjj-1)*ny,ix)=u_elimine{iii,jjj}(1)*speye(ny);end;else;
f=retf(yy,u{iii,jjj},aly);
for ix=1:mx;ff([1:ny]+(iii-1)*ny,[1:ny]+(jjj-1)*ny,ix)=rettoeplitz(f(ix,:),topy);end;
end;end;end;                 % iii jjj
for ix=1:mx;ff(:,:,ix)=inv(ff(:,:,ix));end;
fff=cell(p);[fff{:}]=deal(zeros(nx,nx,ny));
for iii=1:p;for jjj=1:p;fff{iii,jjj}=ff([1:ny]+(iii-1)*ny,[1:ny]+(jjj-1)*ny,:);end;end;
for iii=1:p;for jjj=1:p;fff{iii,jjj}=reshape(retf(xx,reshape(fff{iii,jjj},ny^2,mx),alx),ny,ny,2*nx-1);end;end;
ffff=cell(p);[ffff{:}]=deal(zeros(nx,nx,ny,ny));
%ffff=cell(p);[ffff{:}]=deal(zeros(nx,ny,nx,ny));
for iii=1:p;for jjj=1:p;
% for ii=1:ny;for jj=1:ny;ffff{iii,jjj}(:,:,ii,jj)=rettoeplitz(fff{iii,jjj}(ii,jj,:),topx);end;end;
% %for ii=1:ny;for jj=1:ny;ffff{iii,jjj}(:,ii,:,jj)=rettoeplitz(fff{iii,jjj}(ii,jj,:),topx);end;end;
% ffff{iii,jjj}=permute(ffff{iii,jjj},[1,3,2,4]);
ffff{iii,jjj}=rettoeplitz(fff{iii,jjj},topx,[1,3,2,4]);% modif 3 2014
ffff{iii,jjj}=reshape(ffff{iii,jjj},nx*ny,nx*ny);
end;end;

case 0;   % <*******  Tf(fonction )= ffff *Tf(fonction continue en x et en y)
ffff=cell(1,p);[ffff{:}]=deal(zeros(nx,ny,nx,ny));% u est un vecteur ligne
ff=zeros(ny,ny,mx);fff=zeros(nx,nx,ny);
for iii=1:p;
f=retf(yy,u{iii},aly); % tf en y fonction de x
for ix=1:mx;ff(:,:,ix)=rettoeplitz(f(ix,:),topy);end;
fff=reshape(retf(xx,reshape(ff,ny^2,mx),alx),ny,ny,2*nx-1);
% for ii=1:ny;for jj=1:ny;ffff{iii}(:,:,ii,jj)=rettoeplitz(fff(ii,jj,:),topx);end;end;
% ffff{iii}=permute(ffff{iii},[1,3,2,4]);
ffff{iii}=rettoeplitz(fff,topx,[1,3,2,4]);% Modif 3 2014

%for ii=1:ny;for jj=1:ny;ffff{iii}(:,ii,:,jj)=rettoeplitz(fff(ii,jj,:),topx);end;end;




ffff{iii}=reshape(ffff{iii},nx*ny,nx*ny);
end;
end;     % <*******
ffff=retsparse(ffff);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [a,b]=retinv(a,b);
%inversion a+ib
a=retio(a,-2);b=retio(b,-2);
am=max(a(:))-min(a(:));bm=max(b(:))-min(b(:));
if am>bm;
if bm<10*eps*am;b=spalloc(size(b,1),size(b,2),0);a=retio(inv(a),1);
else
aa=b*inv(a);
a=inv(a+aa*b);b=-a*aa;
a=retio(a,1);b=retio(b,1);
end;

else
if am<10*eps*bm;a=spalloc(size(a,1),size(a,2),0);b=-inv(retio(b));
else 
bb=-a*inv(b);
b=inv(-b+bb*a);a=b*bb;
a=retio(a,1);b=retio(b,1);
end;
end;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [a,ai]=retblkinv(a,decoupe,N,M,n_bloc,ai);
% b=retblkinv(a,decoupe,N,M,n_bloc);
% M*a*N  blocs de longueur decoupe
% M=inv(N);
% b=inv(a);
% si n_bloc est precisé on ne calcule que la contribution du bloc numero n_bloc 
% a=inv(a);
% si ai est precisé inversion economique en mémoire

if nargin<2;a=inv(a);return;end;
if nargin<5;n_bloc=0;end;

n=length(decoupe);
n2=cumsum(decoupe);n1=1+[0,n2(1:n-1)];
if nargin<6;% version normale                  ****************
	if n_bloc==0;
	b=M*a*N;
	a=0;
	for n_bloc=1:n;if n2(n_bloc)>=n1(n_bloc);a=a+N(:,n1(n_bloc):n2(n_bloc))*inv(b(n1(n_bloc):n2(n_bloc),n(n_bloc):n2(n_bloc)))*M(n1(n_bloc):n2(n_bloc),:);end;end;
	else;
	if n2(n_bloc)>=n1(n_bloc);a=N(:,n1(n_bloc):n2(n_bloc))*inv(M(n1(n_bloc):n2(n_bloc),:)*a*N(:,n1(n_bloc):n2(n_bloc)))*M(n1(n_bloc):n2(n_bloc),:);end
	end;	
	
else;       % version economique en memoire    ****************
	if n_bloc==0;
	b=M*retio(a)*N;bi=M*retio(ai)*N;
	a=0;ai=0;
	for n_bloc=1:n;% <<<<<<<<<<<<
	[prv,prvi]=retinv(b(n1(n_bloc):n2(n_bloc),n(n_bloc):n2(n_bloc)),bi(n1(n_bloc):n2(n_bloc),n(n_bloc):n2(n_bloc)));
	a=a+N(:,n1(n_bloc):n2(n_bloc))*retio(prv)*M(n1(n_bloc):n2(n_bloc),:);
	ai=ai+N(:,n1(n_bloc):n2(n_bloc))*retio(prvi)*M(n1(n_bloc):n2(n_bloc),:);
	end;          % <<<<<<<<<<<<
	else;
	[prv,prvi]=retinv(M(n1(n_bloc):n2(n_bloc),:)*retio(a)*N(:,n1(n_bloc):n2(n_bloc)),M(n1(n_bloc):n2(n_bloc),:)*retio(ai)*N(:,n1(n_bloc):n2(n_bloc)));
	a=N(:,n1(n_bloc):n2(n_bloc))*retio(prv)*M(n1(n_bloc):n2(n_bloc),:);
	ai=N(:,n1(n_bloc):n2(n_bloc))*retio(prvi)*M(n1(n_bloc):n2(n_bloc),:);
	end;
	a=retio(a,1);ai=retio(ai,1);
end;        % version economique en memoire ?  ****************

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	function X=pad_zeros(XX,Nx,Ny,Numx,Numy);% dans le cas ou pas assez de points en Nx et Ny
	XX=retio(XX);
    X=zeros(size(Numx));
	num=Numx+(Numy-1)*Nx;
	f=find((Numx<=Nx) & (Numx>=1) & (Numy<=Ny) & (Numy>=1));
	X(f)=XX(num(f));
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [a,d,w]=retcouche_Aurelien(init,w,c);% Cercles Aureliens 
	[p,pp,q,qq,d,ez,hz,fEt,fHt,fDn,fBn,W,pas]=deal([]);
c=real(c);
beta=init{2};nx=init{6};ny=init{7};n1=nx*ny;si=init{8};cao=init{10};pas=init{end}.d;sog=init{end}.sog;
if ~isempty(cao); % changement de coordonnees
kx=retmat(apodisecao(cao{1})*retdiag(beta(1,1:nx)),ny);
ky=retmat(apodisecao(cao{2})*retdiag(beta(2,1:nx:n1)),-nx);
else
kx=retdiag(beta(1,:));ky=retdiag(beta(2,:));
end;

% W={x,y,{Epsilon,Mu,  ZZeta,ZZetac,KKsi,KKsic,Zeta,Ksi}};
% w={epsilon,eepsilon,mu,mmu,W};

[epsilon,eepsilon,mu,mmu,W]=deal(w{:});
%[S,Sc,C,Cc]=deal(W{3}{3:6});
[S,Sc,C,Cc, prv,prv, S_Sc,C_Cc,S_Cc,C_Sc]=deal(W{3}{3:12});

[Nx,Ny]=size(epsilon);
Numx=repmat(rettoeplitz(1:2*nx-1),ny,ny);
topy=rettoeplitz(1:2*ny-1);Numy=topy(repmat(ceil((1:nx*ny)/nx).',1,nx*ny)+ny*repmat((ceil((1:nx*ny)/nx)-1),nx*ny,1));
offsetx=1+Nx/2-nx;offsety=1+Ny/2-ny;
Numx=Numx+offsetx;
Numy=Numy+offsety;
% Num=zeros(size(aa));Num(:)=1:numel(Num);Num=pad_zeros(Num,Nx,Ny,Numx,Numy);
% epsilon=epsilon(Num);eepsilon=eepsilon(Num);
% Num=zeros(size(aa));Num(:)=1:numel(Num);Num=pad_zeros(Num,Nx,Ny,Numx,Numy);
epsilon=pad_zeros(epsilon,Nx,Ny,Numx,Numy);eepsilon=pad_zeros(eepsilon,Nx,Ny,Numx,Numy);
if numel(mu)>1;mu=pad_zeros(mu,Nx,Ny,Numx,Numy);mmu=pad_zeros(mmu,Nx,Ny,Numx,Numy);inv_mmu=inv(mmu);
else;prv=speye(size(epsilon));mmu=(1/mu)*prv;inv_mmu=mu*prv;mu=mu*prv;end;

inv_eepsilon=inv(eepsilon);

%if numel(mu)==1;inv_mmu=inv(mmu)*speye(size(epsilon));mu=mu*speye(size(epsilon));end;


S=pad_zeros(S,Nx,Ny,Numx,Numy);Sc=pad_zeros(Sc,Nx,Ny,Numx,Numy);C=pad_zeros(C,Nx,Ny,Numx,Numy);Cc=pad_zeros(Cc,Nx,Ny,Numx,Numy);

S_Sc=pad_zeros(S_Sc,Nx,Ny,Numx,Numy);C_Cc=pad_zeros(C_Cc,Nx,Ny,Numx,Numy);S_Cc=pad_zeros(S_Cc,Nx,Ny,Numx,Numy);C_Sc=pad_zeros(C_Sc,Nx,Ny,Numx,Numy);


delta_mu=inv_mmu-mu;delta_epsilon=inv_eepsilon-epsilon;

if c~=-2;% inutile pour retbrillouin *********************
if isempty(si);% non utilisation des proprietes de symetrie
ez=[-1i*eepsilon*ky,1i*eepsilon*kx];
hz=[-1i*mmu*ky,1i*mmu*kx];

% a1=[[(kx*eepsilon*ky               +.5*(delta_mu*S_Cc+S_Cc*delta_mu) ) ,...
%     (-kx*eepsilon*kx               +.5*(delta_mu*S_Sc+S_Sc*delta_mu)+mu)];...
%     [(ky*eepsilon*ky               -.5*(delta_mu*C_Cc+C_Cc*delta_mu)-mu)  ,...
%     (-ky*eepsilon*kx               -.5*(delta_mu*C_Sc+S_Cc*delta_mu))]];
% 
% a2=[[(kx*mmu*ky               + .5*(delta_epsilon*S_Cc+S_Cc*delta_epsilon) ) ,...
%     (-kx*mmu*kx               +.5*(delta_epsilon*S_Sc+S_Sc*delta_epsilon)+epsilon)];...
%     [(ky*mmu*ky               -.5*(delta_epsilon*C_Cc+C_Cc*delta_epsilon)-epsilon)  ,...
%     (-ky*mmu*kx               -.5*(delta_epsilon*C_Sc+S_Cc*delta_epsilon) )]];
% 

% 
% a1=[[(kx*eepsilon*ky               +delta_mu*S_Cc ) ,...
%     (-kx*eepsilon*kx               +delta_mu*S_Sc+mu)];...
%     [(ky*eepsilon*ky               -delta_mu*C_Cc-mu)  ,...
%     (-ky*eepsilon*kx               -delta_mu*C_Sc)]];
% 
% a2=[[(kx*mmu*ky               + delta_epsilon*S_Cc ) ,...
%     (-kx*mmu*kx               +delta_epsilon*S_Sc+epsilon)];...
%     [(ky*mmu*ky               -delta_epsilon*C_Cc-epsilon)  ,...
%     (-ky*mmu*kx               -delta_epsilon*C_Sc )]];



a1=[[(kx*eepsilon*ky               +delta_mu*S_Cc ) ,...
    (-kx*eepsilon*kx               +inv_mmu*S_Sc+mu*C_Cc)];...
    [(ky*eepsilon*ky               -inv_mmu*C_Cc-mu*S_Sc)  ,...
    (-ky*eepsilon*kx               -delta_mu*C_Sc)]];

a2=[[(kx*mmu*ky               + delta_epsilon*S_Cc ) ,...
    (-kx*mmu*kx               +inv_eepsilon*S_Sc+epsilon*C_Cc)];...
    [(ky*mmu*ky               -inv_eepsilon*C_Cc-epsilon*S_Sc)  ,...
    (-ky*mmu*kx               -delta_epsilon*C_Sc )]];








else;  % utilisation des proprietes de symetrie

%     a1=(si{2}(:,1:n1)*kx)*eepsilon*(ky*si{3}(1:n1,:)-kx*si{3}(n1+1:2*n1,:))+(si{2}(:,n1+1:2*n1)*ky)*eepsilon*(ky*si{3}(1:n1,:)-kx*si{3}(n1+1:2*n1,:))...                
%     -si{2}(:,n1+1:2*n1)*mu*si{3}(1:n1,:)+si{2}(:,1:n1)*mu*si{3}(n1+1:2*n1,:)...
%     +.5*(   (si{2}(:,1:n1)*S_Cc-si{2}(:,n1+1:2*n1)*C_Cc)*delta_mu*(si{3}(1:n1,:))+(si{2}(:,1:n1)*S_Sc-si{2}(:,n1+1:2*n1)*C_Sc)*delta_mu*(si{3}(n1+1:2*n1,:))...
% +(si{2}(:,1:n1))*delta_mu*(C_Sc*si{3}(1:n1,:)+S_Sc*si{3}(n1+1:2*n1,:))-(si{2}(:,n1+1:2*n1))*delta_mu*(C_Cc*si{3}(1:n1,:)+S_Cc*si{3}(n1+1:2*n1,:))  );
% 
% 
%     a2=(si{4}(:,1:n1)*kx)*mmu*(ky*si{1}(1:n1,:)-kx*si{1}(n1+1:2*n1,:))+(si{4}(:,n1+1:2*n1)*ky)*mmu*(ky*si{1}(1:n1,:)-kx*si{1}(n1+1:2*n1,:))...                
%     -si{4}(:,n1+1:2*n1)*epsilon*si{1}(1:n1,:)+si{4}(:,1:n1)*epsilon*si{1}(n1+1:2*n1,:)...
%     +.5*(   (si{4}(:,1:n1)*S_Cc-si{4}(:,n1+1:2*n1)*C_Cc)*delta_epsilon*(si{1}(1:n1,:))+(si{4}(:,1:n1)*S_Sc-si{4}(:,n1+1:2*n1)*C_Sc)*delta_epsilon*(si{1}(n1+1:2*n1,:))...
% +(si{4}(:,1:n1))*delta_epsilon*(C_Sc*si{1}(1:n1,:)+S_Sc*si{1}(n1+1:2*n1,:))-(si{4}(:,n1+1:2*n1))*delta_epsilon*(C_Cc*si{1}(1:n1,:)+S_Cc*si{1}(n1+1:2*n1,:)) );


    a1=(si{2}(:,1:n1)*kx)*eepsilon*(ky*si{3}(1:n1,:)-kx*si{3}(n1+1:2*n1,:))+(si{2}(:,n1+1:2*n1)*ky)*eepsilon*(ky*si{3}(1:n1,:)-kx*si{3}(n1+1:2*n1,:))...                
    -si{2}(:,n1+1:2*n1)*mu*si{3}(1:n1,:)+si{2}(:,1:n1)*mu*si{3}(n1+1:2*n1,:)...
+(si{2}(:,1:n1))*delta_mu*(C_Sc*si{3}(1:n1,:)+S_Sc*si{3}(n1+1:2*n1,:))-(si{2}(:,n1+1:2*n1))*delta_mu*(C_Cc*si{3}(1:n1,:)+S_Cc*si{3}(n1+1:2*n1,:))  ;


    a2=(si{4}(:,1:n1)*kx)*mmu*(ky*si{1}(1:n1,:)-kx*si{1}(n1+1:2*n1,:))+(si{4}(:,n1+1:2*n1)*ky)*mmu*(ky*si{1}(1:n1,:)-kx*si{1}(n1+1:2*n1,:))...                
    -si{4}(:,n1+1:2*n1)*epsilon*si{1}(1:n1,:)+si{4}(:,1:n1)*epsilon*si{1}(n1+1:2*n1,:)...
+(si{4}(:,1:n1))*delta_epsilon*(C_Sc*si{1}(1:n1,:)+S_Sc*si{1}(n1+1:2*n1,:))-(si{4}(:,n1+1:2*n1))*delta_epsilon*(C_Cc*si{1}(1:n1,:)+S_Cc*si{1}(n1+1:2*n1,:)) ;



ez=-1i*eepsilon*ky*si{3}(1:n1,:)+1i*eepsilon*kx*si{3}(n1+1:2*n1,:);
hz=-1i*mmu*ky*si{1}(1:n1,:)+1i*mmu*kx*si{1}(n1+1:2*n1,:);
end;



a1=1i*a1;a2=-1i*a2;ez=1i*ez;hz=-1i*hz;% declonage car retcouche n'est pas cloné !

end;% **************************************

moharam=0;sens=0;
if c~=-2;[p,d]=retc3(a1,a2,c,moharam,sens);if c~=-3;q=a1*p;else;q=[];end;else p=[];q=[];d=[];end;  % diagonalisation
if c<0;pp=[];qq=[];else pp=retpinv(p);qq=retpinv(q);end;

if c==1;% calcul precis des champs
%  fDn=[inv_eepsilon*pad_zeros(W{3}{6},Nx,Ny,Numx,Numy),inv_eepsilon* pad_zeros(W{3}{4},Nx,Ny,Numx,Numy)];  
%  fEt=[-pad_zeros(W{3}{3},Nx,Ny,Numx,Numy),pad_zeros(W{3}{5},Nx,Ny,Numx,Numy)];
 %W={x,y,{Epsilon,Mu,  ZZeta,ZZetac,KKsi,KKsic,Zeta,Ksi}};

fDn=[inv_eepsilon*Cc,inv_eepsilon*Sc];  
fEt=[-S,C];
if ~isempty(si);
fDn2=[inv_eepsilon*C,inv_eepsilon*S];  
fEt2=[-Sc,Cc];
fDn={fDn,fDn2};
fEt={fEt,fEt2};
end 
 
 
if numel(mu)>1;
fBn=[inv_mmu*Cc,inv_mmu*Sc];  
fHt=[-S,C];
if ~isempty(si);
fBn2=[inv_mmu*C,inv_mmu*S];  
fHt2=[-Sc,Cc];
fBn={fBn,fBn2};
fHt={fHt,fHt2};
end 
 else;fBn=[];fHt=[];end;% mu constant
% W={x,y,{Epsilon,Mu,  ZZeta,ZZetac,KKsi,KKsic,Zeta,Ksi}};

else;
[fDn,fEt,fBn,fHt]=deal({});
end;   
a={p,pp,q,qq,d,ez,hz,fEt,fHt,fDn,fBn,W,pas,init{2},init{11},struct('dim',2,'genre',2,'type',8,'sog',sog)};
%if (nargout==3)& nw;
if (nargout==3);
%w={x/pas(1),y/pas(2),w{3},eepz,mmuz,muy,mux,epy,epx,hx,ez,hz,fex,fhx,fey,fhy};% attention si transformee de granet on remet dans w ev sortie les x y physiques !
                                                                       % Le w dans a a les coordonnees numeriques il peut etre utilise par retchamp
inv_mu=inv(mu);inv_epsilon=inv(epsilon);
%w={eepsilon,mmu,inv_epsilon,inv_mu,S,Sc,C,Cc};

w={eepsilon,mmu,inv_epsilon,inv_mu,S_Sc,C_Cc,S_Cc,C_Sc};

% % muy=ret2li(x,y,u(:,:,2),nx,ny,1,0,li(3));muy=retio(muy,io);  % muy
% % mux=ret2li(x,y,u(:,:,1),nx,ny,0,1,li(4));mux=retio(mux,io);  % mux
% % epy=ret2li(x,y,u(:,:,5),nx,ny,1,0,li(5));epy=retio(epy,io);  % epy
% % epx=ret2li(x,y,u(:,:,4),nx,ny,0,1,li(6));epx=retio(epx,io);  % epx
% % eepz=ret2li(x,y,1./u(:,:,6),nx,ny,0,0,li(7));eepz=retio(eepz,io);
% % mmuz=ret2li(x,y,1./u(:,:,3),nx,ny,0,0,li(8));mmuz=retio(mmuz,io);
end;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   diagonalisation: ATTENTION non cloné 
%  si a1 et a2 sont clonés, a1=1i*a1;a2=-1i*a2;% declonage
function [p,d]=retc3(a1,a2,c,moharam,sens);
if nargin<4;moharam=0;end;
if nargin<3;c=1;end;
a=full(a2*a1);clear a1 a2	
%moharam=0;
if moharam>0;% diagonalisation moharam en 2 blocs
n=size(a,1);
if c==-3;p=[];d=[reteig(a(1:moharam,1:moharam));reteig(a(moharam+1:end,moharam+1:end))];
else;
[p1,d1]=reteig(a(1:moharam,1:moharam));[p2,d2]=reteig(a(moharam+1:end,moharam+1:end));d=diag([diag(d1);diag(d2)]);
if sens>0;% a12=0
pp=((p2\a(moharam+1:end,1:moharam))*p1)./(repmat(diag(d1).',n-moharam,1)-repmat(diag(d2),1,moharam));p=[p1,zeros(moharam,n-moharam);p2*pp,p2];
else;    % a21=0
pp=((p1\a(1:moharam,moharam+1:end))*p2)./(repmat(diag(d2).',moharam,1)-repmat(diag(d1),1,n-moharam));p=[p1,p1*pp;zeros(n-moharam,moharam),p2];
end;
%retcompare(a*p,p*d)
d=diag(d);
end;

else; % diagonalisation complete
if c==-3;d=reteig(a);p=[];else;[p,d]=reteig(a);d=diag(d);end;aa=a;
end;
d=retsqrt(d,1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [a1,a2,k,aa2]=retc2(mu,mmu,ep,beta,cao,a1,aa2);%n=length(beta);
%calcul de la matrice m=[0,a1;a2,0] associee au passage dans un milieu periodique
%invariant en y
%mu :coefficients de fourier de mu (2*n-1 valeurs le coefficient d'ordre 0 etant range en mu(n))
%mmu:idem pour 1/mu   ep:idem pour epsilon
%beta:n constantes de propagation en x du developpement de rayleigh

if length(cao)>0;

%k=apodisecao(cao,0)*diag(beta);% sans apodisation 
k=apodisecao(cao)*diag(beta);% avec apodisation 
else k=diag(beta);end;
if nargin<7;aa2=[];end;if isempty(aa2);
a1=inv(rettoeplitz(mmu));aa2=inv(rettoeplitz(mu));
end;
a2=k*aa2*k-rettoeplitz(ep);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function cao=apodisecao(cao,parm);if nargin<2;parm=1;end;
%parm=0;% test
if iscell(cao);% transformée de Granet
granet=1;
cao_granet=cao{2};cao=cao{1};


else;
granet=0;
end;
if ~isempty(cao);% transformee complexe
n=ceil(length(cao)/2);
if parm==1;%  apodisation    
apod=7+n/1500;if apod>7.5;apod=1;end;
cao(n)=cao(n)-1;
if n<100;pprv=retchamp([apod,n]);
else;
pprv=sum(abs(rettoeplitz(cao)).^2,2)-sum(abs(cao).^2);for ii=1:5;pprv=(pprv([1,1:end-1])+pprv([2:end,end]))/2;end; % lissage
pprv=10.^(-20*sqrt(abs(pprv)));
end;
    if granet;% transformee de Granet
%     cao=conv(cao_granet,cao);cao=cao(n:3*n-2);    
%     cao=diag(pprv)*rettoeplitz(cao)+rettoeplitz(cao_granet);   
    cao=rettoeplitz(cao_granet);% modif 28 8 2011   
    else;
    cao=diag(pprv)*rettoeplitz(cao)+eye(n);
    end;
else;% pas apodisation
    if granet;% transformee de Granet
    cao=conv(cao_granet,cao);cao=cao(n:3*n-2);
    end
    cao=rettoeplitz(cao);
end;% apodisation ?

end;% transformee complexe ?

if granet;% transformee de Granet
if isempty(cao);
cao=rettoeplitz(cao_granet);
% <------------- modif 8 2011
%     else;
%     %cao=rettoeplitz(cao_granet)*cao;    
%     cao=cao*rettoeplitz(cao_granet);    
end;    
end;    

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [fmu,fmmu,fep]=retc1(x,ep,beta);n=length(beta);alpha=(-n+1:n-1)*(2*pi);
fmu=retf(x,ep(2,:),alpha);fmmu=retf(x,ep(3,:),alpha);fep=retf(x,ep(1,:),alpha);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [a,dd]=retincline(init,c,u);
% 1D
% function [a,dd]=retincline(init,X,ep,teta);
%  milieu periodique invariant 
% dans une direction inclinee d'un angle teta(degres) sur y
% ep valeurs de [epsilon;mux;1/muy]
% X={x,wx} points et poids de Gauss 
% maillage de texture: u={{x,wx},ep,teta,struct('type',5)}
% Attention si l'inclinaison est non nulle pas de changement de coordonnees !
% EXEMPLE reseau sinusoidal 
%   mm=40;init=retinit(d,[-mm,mm],beta,sym);
%   dnn=1.e-2;inclinaison=10;
% 	ntr=max(2*mm+1,5);[x,w]=retgauss(0,d,ntr,-4);nn=nn0+dnn*cos(2*pi*x/d);% echantillonage sinusoide
% 	ep=retep(nn,pol,k0);x={x,w};% ep valeurs de [epsilon;mux;1/muy]
% 	u={x,ep,inclinaison*pi/180,struct('type',5)};
% 
% 2D
% function [a,dd]=retincline(init,X,Y,ep,teta,phi);
% X={x,wx} ,Y={y,wy} points et poids de Gauss 
% maillage de texture: u={{x,wx},{y,wy},ep,teta,phi,struct('type',5)}

io=imag(c)~=0;c=real(c);parm=u{end};
if init{end}.dim==1;   % <*****************  1D
x=u{1}{1};ep=u{1}{2};teta=parm.teta;
% [x,ep,teta,parm]=deal(varargin{1:4});
if isfield(parm,'Maystre');Maystre=parm.Maystre;else Maystre=[];end;
if isfield(parm,'fac');ffac=[1;parm.fac(:)];else ffac=[1,1];end;if length(ffac)==1;ffac=[1,1];end;
n=init{1};d=init{end}.d;beta=init{2};si=init{3};cao=init{5};sog=init{end}.sog;n=size(beta,2);
if iscell(x);x{1}=x{1}/d;x{2}=x{2}/d;else x=x/d;end;
if (teta==0)&(size(ep,1)==3);% <.........pas inclinaison et isotrope
[fmu,fmmu,fep]=retc1(x,ep,beta);
[a1,a2,k,aa2]=retc2(fmu,fmmu,fep,beta,cao);
if ~isempty(si);
aa1=a1*si{1};a1=si{2}*aa1;a2=si{2}*a2*si{1};
end    
[p,dd]=retc3(a1,a2,c);q=a1*p;
if c>=0;pp=retpinv(p);qq=retpinv(q);q=a1*p;else q=[];pp=[];qq=[];end

if c==1;
a3=-inv(rettoeplitz(fmu))*k; % matrice a3 permettant le calcul de hy=i/mu*du/dx
if ~isempty(si);a3=a3*si{1};a1=aa1;end
else
a3=[];if ~isempty(si);a1=aa1;end
end;
if iscell(x);xx=[x{1}(2:end),x{1}(1)+1];xx=(x{1}+xx)/2;ww={xx*d,ep};else ww={x*d,ep};end;
a={p,pp,q,qq,dd,a3,a1,ww,d,beta,struct('dim',1,'genre',2,'type',1,'sog',sog)};
	
else  % <...........inclinaison
if length(cao)>0;k=apodisecao(cao)*diag(beta);else k=diag(beta);end;% avec apodisation
if size(ep,1)<4;ep=[ep;0*ep(1,:)];end;% completé par des 0 pour milieux isotropes
ct=cos(teta);st=sin(teta);tt=st/ct;
nn=length(beta);alpha=(-nn+1:nn-1)*(2*pi);
% prv=(st^2)*ep(2,:)+(ct^2)./ep(3,:);
% prv11=retf(x,ep(2,:)./(ep(3,:).*prv),alpha);
% prv12=retf(x,ep(2,:)./prv,alpha);
% prv22=retf(x,1./prv,alpha);

prv=(st^2)*ep(2,:)+(ct^2)./ep(3,:)-2*st*ct*ep(4,:);
prv11=retf(x,ct./prv,alpha);
prv12=retf(x,(st*ep(2,:)-ct*ep(4,:))./prv,alpha);
prv22=retf(x,ct*(ep(2,:)./ep(3,:)-ep(4,:).^2)./prv,alpha);

fep=retf(x,ep(1,:),alpha);prv0=rettoeplitz(fep);

ddd=zeros(length(ffac)-1,2*n);I=speye(nn);O=sparse(nn,nn);
for ifac=1:length(ffac)-1;fac=ffac(ifac);
m=[[rettoeplitz(prv11/fac),rettoeplitz(prv12)];[-rettoeplitz(prv12),rettoeplitz(prv22*fac)]];
%m=(1/ct)*(m\(eye(2*nn)+st*[-m(:,nn+1:2*nn),m(:,1:nn)]));
m=(1/ct)*(inv(m)-st*[[O,-I];[I,O]])*[[O,I];[-i*k,O]];
if ifac==1;
%mm=[[O,st*k];[st*I,O]]+ct*m;mm(:,1:nn)=i*mm(:,1:nn);% declonage
a3=[ct*m(nn+1:2*nn,1:nn),ct*m(nn+1:2*nn,nn+1:2*nn)+st*I];       % matrice a3 permettant le calcul de Ht
a1=[ct*m(1:nn,1:nn)+(i*st)*k,ct*m(1:nn,nn+1:2*nn)];  % matrice a1 permettant le calcul de Bn
end;
m=[[m(1:nn,1:nn)+(i*st/ct)*k,m(1:nn,nn+1:2*nn)];[-prv0*fac+i*k*m(nn+1:2*nn,1:nn),(i*k)*m(nn+1:2*nn,nn+1:2*nn)+(i*st/ct)*k]];
if ifac==1;[p,vp]=reteig(m);p=p*retdiag(1./sqrt(sum(abs(p).^2)));dd=diag(vp);pp=retpinv(p);ddd(1,:)=dd.';    % diagonalisation
else;dddd=pp*m*p;ddd0=eig(dddd);dddd=diag(dddd);ii=retassocie(ddd0,dddd,i);ddd(ifac,:)=ddd0(ii).';
end	
end;
if ifac>1;dddd=interp1(imag(ffac(1:end-1)),ddd,imag(ffac(end)),'pchip').';else;dddd=dd;end;

%figure;hold on;plot(dd,'or');plot(ddd,'*c');plot(dddd,'ob');plot(real([ddd;dddd.']),imag([ddd;dddd.']),'-g');grid
io=tri_incline(dd,dddd,teta,n,Maystre);% tri
p=p(:,io);pp=pp(io,:);dd=dd(io);  % classement
if iscell(x);xx=[x{1}(2:end),x{1}(1)+1];xx=(x{1}+xx)/2;ww={xx*d,ep};else;ww={x*d,ep};end;
a={p,pp,tt,n,dd,a3,a1,ww,repmat(beta,1,2),struct('dim',1,'genre',2,'type',5,'sog',sog,'d',d)};
end;% <...........inclinaison	? ?

else;                 % <*****************  2 D
x=u{1}{1};y=u{1}{2};ep=u{1}{3};teta=parm.teta;phi=parm.phi;
%[x,y,ep,teta,phi,parm]=deal(varargin{1:6});
if isfield(parm,'Maystre');Maystre=parm.Maystre;else Maystre=[];end;
if isfield(parm,'fac');ffac=[1;parm.fac(:)];else;ffac=[1,1];end;if length(ffac)==1;ffac=[1,1];end;
if isfield(parm,'sens');sens=parm.sens;else;sens=1;end;% parametre pour ret2li_matrice 1 cx -->  cy ,  -1  cy --> cx
beta=init{2};d=init{end}.d;nx=init{6};ny=init{7};n1=nx*ny;si=init{8};n=init{1};cao=init{10};sog=init{end}.sog;
ct=cos(teta);st=sin(teta);cp=cos(phi);sp=sin(phi);
vect=struct('TT',[cp,sp]*(st/ct),'n',[st*cp,st*sp,ct],'nx',[0,ct,-st*sp]/sqrt(ct^2+st^2*sp^2),'ny',[ct,0,-st*cp]/sqrt(ct^2+st^2*cp^2));

if iscell(x);mx=size(x{1},2);else;mx=size(x,2);end;% gauss
if iscell(y);my=size(y{1},2);else;my=size(y,2);end;% gauss

eep=ep;ep=zeros(mx,my,3,3,2);for ii=1:3;ep(:,:,ii,ii,1)=eep(:,:,ii);ep(:,:,ii,ii,2)=eep(:,:,ii+3);end;
%if length(size(ep))<4;eep=ep;ep=zeros(mx,my,3,3,2);for ii=1:3;ep(:,:,ii,ii,1)=eep(:,:,ii);ep(:,:,ii,ii,2)=eep(:,:,ii+3);end;end;
if size(eep,3)>6; % anisotrope
for ii=1:2;
[ep(:,:,1,2,ii),ep(:,:,2,1,ii)]=deal(eep(:,:,7+3*(ii-1)));%     1 7 9       4   10  12
[ep(:,:,2,3,ii),ep(:,:,3,2,ii)]=deal(eep(:,:,8+3*(ii-1)));% mu= 7 2 8   ep= 10  5  11
[ep(:,:,1,3,ii),ep(:,:,3,1,ii)]=deal(eep(:,:,9+3*(ii-1)));%     8 9 3       12  11  6
end;
end;
EP=cell(3,3,2);for ii=1:3;for jj=1:3;for kk=1:2;EP{ii,jj,kk}=ep(:,:,ii,jj,kk);end;end;end;
V={1,0,0;0,1,0;vect.n(1),vect.n(2),vect.n(3)};U={1,0,0;0,1,0;-vect.n(1)/vect.n(3),-vect.n(2)/vect.n(3),1/vect.n(3)};
ke=retprod_cell(V(3,:),EP(:,:,2),U);kh=retprod_cell(V(3,:),EP(:,:,1),U);


me=cell(3);me1=cell(3);mh=cell(3);mh1=cell(3);
[me1{:}]=deal(1,ep(:,:,2,1,2)*vect.nx(2)+ep(:,:,3,1,2)*vect.nx(3),vect.n(1),...
			  0,ep(:,:,2,2,2)*vect.nx(2)+ep(:,:,3,2,2)*vect.nx(3),vect.n(2),...
			  0,ep(:,:,2,3,2)*vect.nx(2)+ep(:,:,3,3,2)*vect.nx(3),vect.n(3));
[me{:}]=deal(ep(:,:,1,1,2)*vect.ny(1)+ep(:,:,3,1,2)*vect.ny(3),0,vect.n(1),...
		     ep(:,:,1,2,2)*vect.ny(1)+ep(:,:,3,2,2)*vect.ny(3),1,vect.n(2),...
		     ep(:,:,1,3,2)*vect.ny(1)+ep(:,:,3,3,2)*vect.ny(3),0,vect.n(3));
[mh1{:}]=deal(1,ep(:,:,2,1,1)*vect.nx(2)+ep(:,:,3,1,1)*vect.nx(3),vect.n(1),...
			 0,ep(:,:,2,2,1)*vect.nx(2)+ep(:,:,3,2,1)*vect.nx(3),vect.n(2),...
			 0,ep(:,:,2,3,1)*vect.nx(2)+ep(:,:,3,3,1)*vect.nx(3),vect.n(3));
[mh{:}]=deal(ep(:,:,1,1,1)*vect.ny(1)+ep(:,:,3,1,1)*vect.ny(3),0,vect.n(1),...
			 ep(:,:,1,2,1)*vect.ny(1)+ep(:,:,3,2,1)*vect.ny(3),1,vect.n(2),...
			 ep(:,:,1,3,1)*vect.ny(1)+ep(:,:,3,3,1)*vect.ny(3),0,vect.n(3));
		 
if sens==1;  %  < SENS = 1
me=retprod_cell(me1,retinv_cell(me));mh=retprod_cell(mh1,retinv_cell(mh));
if c==1;  % calcul precis des champs
f_champe=cell(1,my);[f_champe{:}]=deal(zeros(3*n1));for ii=1:my;prv=me;for jj=1:numel(prv);prv{jj}=prv{jj}(:,ii);end;f_champe{ii}=ret2li_matrice(x,1,prv,nx,ny,mx,1,sens);end;f_champe=retio(f_champe,io);
f_champh=cell(1,my);[f_champh{:}]=deal(zeros(3*n1));for ii=1:my;prv=mh;for jj=1:numel(prv);prv{jj}=prv{jj}(:,ii);end;f_champh{ii}=ret2li_matrice(x,1,prv,nx,ny,mx,1,sens);end;f_champh=retio(f_champh,io);
else;f_champe=[];f_champh=[];
end;     % calcul precis des champs ?
else;       %  < SENS = -1
me=retprod_cell(me,retinv_cell(me1));mh=retprod_cell(mh,retinv_cell(mh1));
if c==1;  % calcul precis des champs
f_champe=cell(1,mx);[f_champe{:}]=deal(zeros(3*n1));for ii=1:mx;prv=me;for jj=1:numel(prv);prv{jj}=prv{jj}(ii,:);end;f_champe{ii}=ret2li_matrice(1,y,prv,nx,ny,1,my,sens);end;f_champe=retio(f_champe,io);
f_champh=cell(1,mx);[f_champh{:}]=deal(zeros(3*n1));for ii=1:mx;prv=mh;for jj=1:numel(prv);prv{jj}=prv{jj}(ii,:);end;f_champh{ii}=ret2li_matrice(1,y,prv,nx,ny,1,my,sens);end;f_champh=retio(f_champh,io);
else;f_champe=[];f_champh=[];
end;     % calcul precis des champs ?
end;        %  < SENS ?	
I=speye(n1);O=sparse(n1,n1);
if ~isempty(cao); % changement de coordonnees
kx=retmat(apodisecao(cao{1})*retdiag(beta(1,1:nx)),ny);
ky=retmat(apodisecao(cao{2})*retdiag(beta(2,1:nx:n1)),-nx);
else
kx=retdiag(beta(1,:));ky=retdiag(beta(2,:));
end;
me=ret2li_matrice(x,y,me,nx,ny,mx,my,sens);
mh=ret2li_matrice(x,y,mh,nx,ny,mx,my,sens);
ke=ret2li_matrice(x,y,ke,nx,ny,mx,my,0);ke=retprod_cell(ke,V);
kh=ret2li_matrice(x,y,kh,nx,ny,mx,my,0);kh=retprod_cell(kh,V);

ddd=zeros(length(ffac)-1,2*n);
for ifac=1:length(ffac)-1;fac=ffac(ifac);      % <........ifac
if sens==1;  %  < SENS = 1
ne=[(-vect.n(3)*fac)*me{1,3},vect.ny(1)*I,-vect.nx(2)*me{1,2};...
	-vect.n(3)*me{2,3},O,(-vect.nx(2)/fac)*me{2,2};...
	(-fac)*ke{1,3},vect.n(1)*I,vect.n(2)*I];
ne_prv=[me{1,1}*fac+(vect.n(1)*fac)*me{1,3},(vect.n(2)*fac)*me{1,3},((-i*vect.nx(3))*me{1,2}+(i*vect.ny(3))*I)*ky,   i*(vect.nx(3)*me{1,2}-vect.ny(3)*I)*kx;...
	me{2,1}+vect.n(1)*me{2,3},vect.n(2)*me{2,3}-I,-(i*vect.nx(3)/fac)*me{2,2}*ky,(i*vect.nx(3)/fac)*me{2,2}*kx;...
	ke{1,1}*fac,ke{1,2}*fac,(i*vect.n(3))*ky,(-i*vect.n(3))*kx];
try;[l_prv,u_prv,p_prv,q_prv]=lu(sparse(ne));ne=(q_prv/u_prv)*(l_prv\(p_prv*ne_prv));clear l_prv u_prv p_prv q_prv ne_prv
catch;[l_prv,u_prv,p_prv]=lu(sparse(ne));ne=u_prv\(l_prv\(p_prv*ne_prv));clear l_prv u_prv p_prv ne_prv
end;	
% ne=retfullsparse([(-vect.n(3)*fac)*me{1,3},vect.ny(1)*I,-vect.nx(2)*me{1,2};...
% 	-vect.n(3)*me{2,3},O,(-vect.nx(2)/fac)*me{2,2};...
% 	(-fac)*ke{1,3},vect.n(1)*I,vect.n(2)*I])\...
% 	[me{1,1}*fac+(vect.n(1)*fac)*me{1,3},(vect.n(2)*fac)*me{1,3},((-i*vect.nx(3))*me{1,2}+(i*vect.ny(3))*I)*ky,   i*(vect.nx(3)*me{1,2}-vect.ny(3)*I)*kx;...
% 	me{2,1}+vect.n(1)*me{2,3},vect.n(2)*me{2,3}-I,-(i*vect.nx(3)/fac)*me{2,2}*ky,(i*vect.nx(3)/fac)*me{2,2}*kx;...
% 	ke{1,1}*fac,ke{1,2}*fac,(i*vect.n(3))*ky,(-i*vect.n(3))*kx];

nh=[(-vect.n(3)*fac)*mh{1,3}, vect.ny(1)*I,(-vect.nx(2))*mh{1,2};...
	(-vect.n(3))*mh{2,3},O,(-vect.nx(2)/fac)*mh{2,2};...
	(-fac)*kh{1,3},vect.n(1)*I,vect.n(2)*I];
nh_prv=[((-i*vect.nx(3))*mh{1,2}+(i*vect.ny(3))*I)*ky,((i*vect.nx(3))*mh{1,2}-(i*vect.ny(3))*I)*kx,fac*mh{1,1}+(vect.n(1)*fac)*mh{1,3},(vect.n(2)*fac)*mh{1,3};...
	 (-i*vect.nx(3)/fac)*mh{2,2}*ky,(i*vect.nx(3)/fac)*mh{2,2}*kx,mh{2,1}+vect.n(1)*mh{2,3},vect.n(2)*mh{2,3}-I ;...
	(i*vect.n(3))*ky,(-i*vect.n(3))*kx, kh{1,1}*fac,kh{1,2}*fac];
try;[l_prv,u_prv,p_prv,q_prv]=lu(sparse(nh));nh=(q_prv/u_prv)*(l_prv\(p_prv*nh_prv));clear l_prv u_prv p_prv q_prv nh_prv
catch;[l_prv,u_prv,p_prv]=lu(sparse(nh));nh=u_prv\(l_prv\(p_prv*nh_prv));clear l_prv u_prv p_prv  nh_prv
end;
% nh=retfullsparse([(-vect.n(3)*fac)*mh{1,3}, vect.ny(1)*I,(-vect.nx(2))*mh{1,2};...
% 	(-vect.n(3))*mh{2,3},O,(-vect.nx(2)/fac)*mh{2,2};...
% 	(-fac)*kh{1,3},vect.n(1)*I,vect.n(2)*I])\...
% 	[((-i*vect.nx(3))*mh{1,2}+(i*vect.ny(3))*I)*ky,((i*vect.nx(3))*mh{1,2}-(i*vect.ny(3))*I)*kx,fac*mh{1,1}+(vect.n(1)*fac)*mh{1,3},(vect.n(2)*fac)*mh{1,3};...
% 	 (-i*vect.nx(3)/fac)*mh{2,2}*ky,(i*vect.nx(3)/fac)*mh{2,2}*kx,mh{2,1}+vect.n(1)*mh{2,3},vect.n(2)*mh{2,3}-I ;...
% 	(i*vect.n(3))*ky,(-i*vect.n(3))*kx, kh{1,1}*fac,kh{1,2}*fac];
if ifac==1;if c==1;  % calcul precis des champs attention si ifac==1 on ne recalcule pas Ze et Zh
Ze=retio([[I,O,O,O];[O,O,-i*vect.nx(3)*ky,i*vect.nx(3)*kx]+vect.nx(2)*ne(2*n1+1:3*n1,:);[vect.n(1)*I,vect.n(2)*I,O,O]+vect.n(3)*ne(1:n1,:)],io);
Zh=retio([[O,O,I,O];[-i*vect.nx(3)*ky,i*vect.nx(3)*kx,O,O]+vect.nx(2)*nh(2*n1+1:3*n1,:);[O,O,vect.n(1)*I,vect.n(2)*I]+vect.n(3)*nh(1:n1,:)],io);
else;Ze=[];Zh=[];
end;end;     % calcul precis des champs ?

else;       %  < SENS = -1
	
ne=[-vect.n(3)*me{1,3},-(vect.ny(1)/fac)*me{1,1},O;...
	-(vect.n(3)*fac)*me{2,3},-vect.ny(1)*me{2,1},vect.nx(2)*I;...
	(-fac)*ke{1,3},vect.n(1)*I,vect.n(2)*I];
ne_prv=[vect.n(1)*me{1,3}-I,vect.n(2)*me{1,3}+me{1,2},(-i*vect.ny(3)/fac)*me{1,1}*ky,(i*vect.ny(3)/fac)*me{1,1}*kx;...
	(vect.n(1)*fac)*me{2,3},(vect.n(2)*fac)*me{2,3}+fac*me{2,2},-((i*vect.ny(3))*me{2,1}-i*vect.nx(3)*I)*ky,((i*vect.ny(3))*me{2,1}-i*vect.nx(3)*I)*kx;...
	ke{1,1}*fac,ke{1,2}*fac,(i*vect.n(3))*ky,(-i*vect.n(3))*kx];


try;[l_prv,u_prv,p_prv,q_prv]=lu(sparse(ne));ne=(q_prv/u_prv)*(l_prv\(p_prv*ne_prv));clear l_prv u_prv p_prv q_prv ne_prv
catch;[l_prv,u_prv,p_prv]=lu(sparse(ne));ne=u_prv\(l_prv\(p_prv*ne_prv));clear l_prv u_prv p_prv ne_prv
end;

% ne=retfullsparse([-vect.n(3)*me{1,3},-(vect.ny(1)/fac)*me{1,1},O;...
% 	-(vect.n(3)*fac)*me{2,3},-vect.ny(1)*me{2,1},vect.nx(2)*I;...
% 	(-fac)*ke{1,3},vect.n(1)*I,vect.n(2)*I])\...
% 	[vect.n(1)*me{1,3}-I,vect.n(2)*me{1,3}+me{1,2},(-i*vect.ny(3)/fac)*me{1,1}*ky,(i*vect.ny(3)/fac)*me{1,1}*kx;...
% 	(vect.n(1)*fac)*me{2,3},(vect.n(2)*fac)*me{2,3}+fac*me{2,2},-((i*vect.ny(3))*me{2,1}-i*vect.nx(3)*I)*ky,((i*vect.ny(3))*me{2,1}-i*vect.nx(3)*I)*kx;...
% 	ke{1,1}*fac,ke{1,2}*fac,(i*vect.n(3))*ky,(-i*vect.n(3))*kx];
nh=[-vect.n(3)*mh{1,3},-(vect.ny(1)/fac)*mh{1,1},O;...
	(-vect.n(3)*fac)*mh{2,3},-vect.ny(1)*mh{2,1},vect.nx(2)*I;...
	(-fac)*kh{1,3},vect.n(1)*I,vect.n(2)*I];
nh_prv=[(-i*vect.ny(3)/fac)*mh{1,1}*ky,(i*vect.ny(3)/fac)*mh{1,1}*kx,vect.n(1)*mh{1,3}-I,vect.n(2)*mh{1,3}+mh{1,2};...
	-((i*vect.ny(3))*mh{2,1}-(i*vect.nx(3))*I)*ky,((i*vect.ny(3))*mh{2,1}-i*vect.nx(3)*I)*kx,(vect.n(1)*fac)*mh{2,3},(vect.n(2)*fac)*mh{2,3}+fac*mh{2,2};...
	(i*vect.n(3))*ky,(-i*vect.n(3))*kx, kh{1,1}*fac,kh{1,2}*fac];
try;[l_prv,u_prv,p_prv,q_prv]=lu(sparse(nh));nh=(q_prv/u_prv)*(l_prv\(p_prv*nh_prv));clear l_prv u_prv p_prv q_prv nh_prv
catch;[l_prv,u_prv,p_prv]=lu(sparse(nh));nh=u_prv\(l_prv\(p_prv*nh_prv));clear l_prv u_prv p_prv nh_prv
end;
% nh=retfullsparse([-vect.n(3)*mh{1,3},-(vect.ny(1)/fac)*mh{1,1},O;...
% 	(-vect.n(3)*fac)*mh{2,3},-vect.ny(1)*mh{2,1},vect.nx(2)*I;...
% 	(-fac)*kh{1,3},vect.n(1)*I,vect.n(2)*I])\...
% 	[(-i*vect.ny(3)/fac)*mh{1,1}*ky,(i*vect.ny(3)/fac)*mh{1,1}*kx,vect.n(1)*mh{1,3}-I,vect.n(2)*mh{1,3}+mh{1,2};...
% 	-((i*vect.ny(3))*mh{2,1}-(i*vect.nx(3))*I)*ky,((i*vect.ny(3))*mh{2,1}-i*vect.nx(3)*I)*kx,(vect.n(1)*fac)*mh{2,3},(vect.n(2)*fac)*mh{2,3}+fac*mh{2,2};...
% 	(i*vect.n(3))*ky,(-i*vect.n(3))*kx, kh{1,1}*fac,kh{1,2}*fac];
if ifac==1;if c==1;  % calcul precis des champs
Ze=retio([[O,O,-i*vect.ny(3)*ky,i*vect.ny(3)*kx]+vect.ny(1)*ne(n1+1:2*n1,:);[O,I,O,O];[vect.n(1)*I,vect.n(2)*I,O,O]+vect.n(3)*ne(1:n1,:)],io);
Zh=retio([[-i*vect.ny(3)*ky,i*vect.ny(3)*kx,O,O]+vect.ny(1)*nh(n1+1:2*n1,:);[O,O,O,I];[O,O,vect.n(1)*I,vect.n(2)*I]+vect.n(3)*nh(1:n1,:)],io);
else;Ze=[];Zh=[];
end;end;     % calcul precis des champs ?

end;        %  < SENS ?	

% m=[nh(2*n1+1:3*n1,:)+i*kx*ne(1:n1,:);-nh(n1+1:2*n1,:)+i*ky*ne(1:n1,:);ne(2*n1+1:3*n1,:)+i*kx*nh(1:n1,:);-ne(n1+1:2*n1,:)+i*ky*nh(1:n1,:)]...
% 	+i*vect.TT(1)*diag(repmat(beta(1,:),1,4))+i*vect.TT(2)*diag(repmat(beta(2,:),1,4));
clear kh ke eep ep EP;
ne=sparse(ne);nh=sparse(nh);
ne=retfullsparse([nh(2*n1+1:3*n1,:)+i*kx*ne(1:n1,:);-nh(n1+1:2*n1,:)+i*ky*ne(1:n1,:);ne(2*n1+1:3*n1,:)+i*kx*nh(1:n1,:);-ne(n1+1:2*n1,:)+i*ky*nh(1:n1,:)]...
	+(i*vect.TT(1))*retmat(kx,4)+(i*vect.TT(2))*retmat(ky,4));
clear nh;
ne(1:2*n1,2*n1+1:4*n1)=i*ne(1:2*n1,2*n1+1:4*n1);% declonage
ne(2*n1+1:4*n1,1:2*n1)=-i*ne(2*n1+1:4*n1,1:2*n1);

if ~isempty(si);% utilisation des proprietes de symetrie
ne=[si{2},sparse(n,2*n1);sparse(n,2*n1),si{4}]*ne*[si{1},sparse(2*n1,n);sparse(2*n1,n),si{3}];
end;	


if ifac==1;	
moharam=(nnz(ne(1:n,1:n))==0)&(nnz(ne(n+1:2*n,n+1:2*n))==0);
if moharam; % <*** ne se separe en 2 matrices a1 et a2
a1=ne(1:n,n+1:2*n);a2=ne(n+1:2*n,1:n);
[p,dd]=retc3(a1,a2,c);
if c~=-3;
q=a1*p*diag(1./dd);
pp=retpinv(p);qq=retpinv(q);
[p,pp]=deal([q,q;-p,p],.5*([qq,-pp;qq,pp]));
else;pp=[];end;% seulement les valeurs propres
dd=[-dd;dd];dddd=dd;
else;	   % <***  diagonalisation
if c~=-3;	
[p,vp]=reteig(ne);p=p*retdiag(1./sqrt(sum(abs(p).^2)));pp=retpinv(p);                       
dd=diag(vp);
else;dd=reteig(ne);p=[];q=[];end;
end;      % <***  ne se separe en 2 matrices ?
ddd(1,:)=dd.';
else;dddd=pp*ne*p;ddd0=eig(dddd);dddd=diag(dddd);ii=retassocie(ddd0,dddd,i);ddd(ifac,:)=ddd0(ii).';
end;
if moharam break;end;
end;                                              % <........ifac

%if ifac>1;dddd=interp1(imag(ffac(1:end-1)),ddd,imag(ffac(end)),'pchip').';else;dddd=dd;end;
if ifac>1;dddd=ddd(end-1,:).'*((ffac(end)-ffac(end-1))/(ffac(end-2)-ffac(end-1)))+ddd(end,:).'*((ffac(end)-ffac(end-2))/(ffac(end-1)-ffac(end-2)));else;dddd=dd;end;
%if ifac>1;figure;hold on;plot(dd,'.k');
%	figure;hold on;plot(dd,'or');plot(ddd,'.k');plot(dddd,'*b');plot(real([ddd;dddd.']),imag([ddd;dddd.']),'-g');grid;stop;
%end;
io=tri_incline(dd,dddd,teta,n,Maystre);if c~=-3;p=p(:,io);pp=pp(io,:);end;dd=dd(io);% tri

if iscell(x);xx=[x{1}(2:end),x{1}(1)+1];x=(x{1}+xx)/2;end;% cas de gauss
if iscell(y);yy=[y{1}(2:end),y{1}(1)+1];y=(y{1}+yy)/2;end;
a={p,pp,vect.TT,n,dd,f_champe,f_champh,vect,init{11},Ze,Zh,{x,y,u{1}{3}},init{end}.d,beta,struct('dim',init{end}.dim,'genre',2,'type',5,'sog',sog,'sens',sens)};


if teta==0% <.........pas inclinaison
else;	% <...........inclinaison
end;    % <...........inclinaison	??

end;                  % <*****************  1 D 2 D ?
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function io=tri_incline(dd,dddd,teta,n,Maystre);dd=dd(:);dddd=dddd(:);
if isempty(Maystre);dddd=real(dddd)*sin(abs(teta)/2+pi/4)-imag(dddd)*cos(abs(teta)/2+pi/4);
else;dddd=real(dddd)*cos(Maystre)-imag(dddd)*sin(Maystre);end;% a verifier
[fp,fm]=retfind(dddd>0);
[prv,io]=sort(dddd(fp));fp=fp(io);[prv,io]=sort(-dddd(fm));fm=fm(io);
io=[fm;fp];
ii=retassocie(dd(io(1:n)),-dd(io(n+1:2*n)),i);io(n+1:2*n)=io(n+ii);     % association par couples
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%                  CYLINDRES POPOV                       %%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [a,d]=retcouche_popov(init,u,c);
N=init{end}.nhanckel;L=init{2};k=init{3};wk=init{4};Psip=retio(init{7}{1});Psim=retio(init{7}{2});cao=init{5};
Pml=init{10};
r=u{1};ep=u{2};
r=[r,0];

if ~isempty(Pml);rr=retelimine([r,Pml{1}(:,1).']);else;rr=r;end;% on ajoute les pml reelles aux discontinuites
cc=10;
if length(rr)==1;x=zeros(1,0);wx=zeros(1,0);else;[x,wx]=retgauss(0,max(rr),cc,ceil(4*max(k)*max(rr)/cc),rr);end;
%if cao(2)~=1;[xx,wxx]=retgauss(cao(1),cao(3),cc,ceil(2*max(k)*cao(3)/cc));x=[x,xx];wx=[wx,wxx];end;% si pml

kx=k.'*x;
JL=besselj(L,kx);JLp1=besselj(L+1,kx);JLm1=besselj(L-1,kx);% ne pas mettre retbessel !

I=speye(N);K=retdiag(k);W=retdiag(wk);

[a2,Hz,epErp,epErm]=popov(I,K,W,Psip,Psim,r,ep(4:6,:),ep(1:3,:),x,wx,JL,JLp1,JLm1,cao,Pml,c);
[a1,Ez,muHrp,muHrm]=popov(I,K,W,Psip,Psim,r,ep(1:3,:),ep(4:6,:),x,wx,JL,JLp1,JLm1,cao,Pml,c);
% full(a1),full(a2)
[p,d]=retc3(a1,a2,c);
if c<0;
pp=[];qq=[];q=[];
else;
q=a1*p;pp=retpinv(p);qq=retpinv(q);
end;

a={p,pp,q,qq,d,{Ez,epErp,epErm,Hz,muHrp,muHrm,Psim,Psip},u,struct('dim',1,'genre',2,'type',6,'sog',init{end}.sog)};
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [A,Ez,Ap,Am]=popov(I,K,W,Psip,Psim,r,mu,ep,x,wx,JL,JLp1,JLm1,cao,Pml,c);
[F,f_inf]=calF(x,r,1./mu(1,:),cao,Pml,1);
Ap=inv(f_inf*I+calFF(JLp1,JLp1,x,wx,K,W,F));
Am=inv(f_inf*I+calFF(JLm1,JLm1,x,wx,K,W,F));

[F,f_inf]=calF(x,r,mu(2,:),cao,Pml,2);
Bp=f_inf*I+calFF(JLp1,JLp1,x,wx,K,W,F);
Bm=f_inf*I+calFF(JLm1,JLm1,x,wx,K,W,F);
Cp=-f_inf*Psim+calFF(JLp1,JLm1,x,wx,K,W,F);
Cm=-f_inf*Psip+calFF(JLm1,JLp1,x,wx,K,W,F);

[F,f_inf]=calF(x,r,ep(3,:),cao,Pml,3);
Ez=inv(f_inf*I+calFF(JL,JL,x,wx,K,W,F));

D=K*Ez*K;if c~=1;Ez=[];end;
A=-.5*i*[[Ap+Bp-D,Ap*Psim+Cp+D];[-Am*Psip-Cm-D,-Am-Bm+D]];
if c~=1;Ap=[];Am=[];end;% pas calcul des champs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [F,f_inf]=calF(x,r,f,cao,Pml,pol);
f_inf=f(1);% valeur externe 
if isempty(x);F=0;return;end;
%if (max(abs(f))==0) & (cao(2)==1);F=0;return;end;
F=f_inf*ones(size(x));

for ii=1:length(r)-1;jj=find((x>r(ii+1)) & (x<r(ii)));F(jj)=f(ii+1);end;
if ~isempty(Pml);
g=retinterp_popov(x,Pml,1)./x;
gg=retinterp_popov(x,Pml,3);

gg_inf=cao(2);

switch pol;
case 1;F=F.*gg./g;   % r
case 2;F=F.*gg./g;   % teta
case 3;F=F.*gg.*g;f_inf=f_inf.*gg_inf^2;   % z
end;
end;                           % << fin pml 

% if cao(2)~=1;%<<<<  pml ( r a ete prolongé jusqu'à r1 dans retu )
% r0=cao(1);r1=cao(3);pml=cao(2);
% rr0=cao(4);ppml=cao(5);
% jj=find(x>r0);F(jj)=f_inf;
% %P=([[3*r0^2,2*r0,1,0];[3*r1^2,2*r1,1,0];[r0^3,r0^2,r0,1];[r1^3,r1^2,r1,1]]\[ppml;pml;rr0;pml*r1]).';
% P=([[20*r0^3,12*r0^2,6*r0,2,0,0];[20*r1^3,12*r1^2,6*r1,2,0,0];[5*r0^4,4*r0^3,3*r0^2,2*r0,1,0];[5*r1^4,4*r1^3,3*r1^2,2*r1,1,0];[r0^5,r0^4,r0^3,r0^2,r0,1];[r1^5,r1^4,r1^3,r1^2,r1,1]]\[0;0;ppml;pml;rr0;pml*r1]).';
% %P=[r1*pml-rr0,rr0*r1*(1-pml)]/(r1-r0);
% PP=polyder(P);g=polyval(P,x(jj))./x(jj);gg=polyval(PP,x(jj));%g_inf=polyval(P,r1)./r1;gg_inf=polyval(PP,r1);
% switch pol;
% case 1;F(jj)=F(jj).*gg./g;   % r
% case 2;F(jj)=F(jj).*gg./g;   % teta
% case 3;F(jj)=F(jj).*gg.*g;f_inf=f_inf.*pml^2;   % z
% end;
% end;%<<<<  pml
% 
%figure;plot(x,real(F),'-k',x,imag(F),'-r');legend('real','imag');title(rettexte(pol,f_inf));

F=F-f_inf;% on retranche la valeur externe
%figure;plot(x,real(F),'.-k',x,imag(F),'.-r');legend('real','imag');title(rettexte(pol,f_inf));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function F=calFF(JL1,JL2,x,wx,K,W,F);
if isempty(x);F=0;return;end;
F=JL1*retdiag(F.*x.*wx)*JL2.';
F=F*K*W;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%                  ELEMENTS FINIS    1 D                 %%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [a_elf,m,me,Maillage]=retelf_couche(init,w,parm,k0,Maillage)
if init{end}.dim==2;[a_elf,m,me,Maillage]=retelf_couche_2D(init,w,parm,k0,Maillage);return;end; % 2D
cao=init{end}.cao;d=init{end}.d;beta0=init{end}.beta0;sym=init{end}.sym;sog=init{end}.sog;if isempty(sym);sym=[0,0];end;
if isempty(cao);cao=[d/2,0];end;xcao=cao(1);lcao=real(cao(2));caoi=imag(cao(2));
if isempty(Maillage);Maillage=retmaillage(parm,w,d,sym,cao);end;  %  maillage
x=Maillage.CoorN;maille=Maillage.Cn;indices=Maillage.NumDom;mh=Maillage.mh;mb=Maillage.mb;msym=Maillage.msym;mbord=Maillage.mbord;
ep=Maillage.ep;

ep(:,1:2)=ep(:,1:2)*k0;ep(:,3)=ep(:,3)/k0;
fac=exp(i*beta0*d);

m=size(x,1);me=size(maille,1);

%%%%%%%%%%%%%%%%%%%%%%
%   ELEMENTS FINIS   %
%%%%%%%%%%%%%%%%%%%%%%
xx=x(abs(maille),:);f=find(maille<0);xx(f,1)=xx(f,1)+d;
xx=reshape(xx,me,3,2);
% surfaces des elements triangulaires et  gradients
%ss=.5*abs((xx(:,2,1)-xx(:,1,1)).*(xx(:,3,2)-xx(:,1,2))-(xx(:,3,1)-xx(:,1,1)).*(xx(:,2,2)-xx(:,1,2)));
ss=(xx(:,2,1)-xx(:,1,1)).*(xx(:,3,2)-xx(:,1,2))-(xx(:,3,1)-xx(:,1,1)).*(xx(:,2,2)-xx(:,1,2));
grad_grad=zeros(me,3,3);
grad_grad(:,1,1)=(xx(:,2,2)-xx(:,3,2))./ss;grad_grad(:,1,2)=(xx(:,3,2)-xx(:,1,2))./ss;grad_grad(:,1,3)=(xx(:,1,2)-xx(:,2,2))./ss;
grad_grad(:,2,1)=(xx(:,3,1)-xx(:,2,1))./ss;grad_grad(:,2,2)=(xx(:,1,1)-xx(:,3,1))./ss;grad_grad(:,2,3)=(xx(:,2,1)-xx(:,1,1))./ss;
grad_grad(:,3,1)=(xx(:,2,1).*xx(:,3,2)-xx(:,3,1).*xx(:,2,2))./ss;grad_grad(:,3,2)=(xx(:,3,1).*xx(:,1,2)-xx(:,1,1).*xx(:,3,2))./ss;
grad_grad(:,3,3)=(xx(:,1,1).*xx(:,2,2)-xx(:,2,1).*xx(:,1,2))./ss;
ss=.5*abs(ss);
% CHANGEMENT DE COORDONNEES

g=1-1/(1+i*caoi);    
A=[(.5-g/8)^2+1/8+g^2/128;g/16-.5;(1+g*(1-g/4))/8;-g/16;(g^2)/128];
xxx=xx;xxx(:,:,1)=mod(xxx(:,:,1)-xcao+d/2,d)-d/2;
offset=mean(xx(:,:,1)-xxx(:,:,1),2);
fcao=find(all(abs(xxx(:,:,1))<=(.5+100*eps)*lcao,2));

offset=offset(fcao);
% calcul des ft des triangles et de leurs derivees
AA=ss;BB=zeros(me,1);CC=zeros(me,1);
if lcao>0;
xxx=xxx(fcao,:,:);
u=(0:4)/lcao;v=zeros(1,5);A1=A(2:5).*(1:4).'/lcao;
[tf_faces,prv]=rettfsimplex(xxx(:,:,1).',xxx(:,:,2).',u,v);
du_tf_faces=prv(:,:,2:5,1);dv_tf_faces=prv(:,:,2:5,2);

AA(fcao)=squeeze(real(tf_faces(:,1,:)))*A+.5*squeeze(real(du_tf_faces(:,1,:)))*A1;
BB(fcao)=.5*squeeze(real(dv_tf_faces(:,1,:)))*A1;
CC(fcao)=squeeze(pi*imag(tf_faces(:,1,2:5)))*A1;
AA(fcao)=AA(fcao)+offset.*CC(fcao);
end;
% 
clear tf_faces du_tf_faces dv_tf_faces xx xxx;

%  CONSTRUCTION DE LA MATRICE CREUSE DES ELEMENTS FINIS
aa=zeros(me,3,3);

integrales=[2,1,1;1,2,1;1,1,2]/12;
for ii=1:3;for jj=1:3;
aa(:,ii,jj)=-ep(indices,1).*ss*integrales(ii,jj)...
+(AA.*grad_grad(:,1,ii).*grad_grad(:,1,jj)+BB.*grad_grad(:,2,ii).*grad_grad(:,1,jj)+CC.*grad_grad(:,3,ii).*grad_grad(:,1,jj))./ep(indices,2)...
+ss.*grad_grad(:,2,ii).*ep(indices,3).*grad_grad(:,2,jj);
end;end;
aa=permute(aa,[2,3,1]);
% pseudo periodicite
periodicite=find(~all(maille>0,2)).';
for ii=periodicite;
k=maille(ii,:);f=find(k<0);
aaa=ones(1,3);aaa(f)=fac;
aa(:,:,ii)=aa(:,:,ii).*((1./aaa.')*aaa);
end;

clear  ss AA BB CC fcao offset ;

% construction de la matrice a
k1=zeros(3,3,me);k2=zeros(3,3,me);
for ii=1:3;prv=repmat(abs(maille(:,ii)).',3,1);
k1(ii,:,:)=reshape(prv,1,3,me);
k2(:,ii,:)=reshape(prv,3,1,me);;
end;
a=sparse(k1(:),k2(:),aa(:),m,m,9*me);
clear aa ;

% construction des matrices Hx Hy Bx qui permettent le calcul des champs
k1=repmat(1:me,1,3);
for ii=periodicite;
f=find(maille(ii,:)<0);
grad_grad(ii,:,f)=grad_grad(ii,:,f)*fac;
end;
% Hx= -i dE/dy 1/mux
g=-i*grad_grad(:,2,:).*reshape(repmat(ep(indices,3),3,1),me,1,3);
Hx=sparse(k1(:),abs(maille(:)),g(:),me,m,9*me);
% Hy= i dE/dx   1/muy
g=i*grad_grad(:,1,:)./reshape(repmat(ep(indices,2),3,1),me,1,3);
Hy=sparse(k1(:),abs(maille(:)),g(:),me,m,9*me);
% Bx= -i dE/dy = mux * Hx
g=-i*grad_grad(:,2,:);
Bx=sparse(k1(:),abs(maille(:)),g(:),me,m,9*me);

clear g k1 k2 prv  ;



%%%%%%%%w%%%%%%%%%%%%%%%%%%%
%  ELIMINATION  :aa*E=H   %
% %%%%%%%%%%%%%%%%%%%%%%%%%%%
% 	% test de symetrie
% 	aaa=a(1:2*msym,1:2*msym);
% 	max(max(abs(aaa(1:msym,1:msym)-aaa(msym+1:2*msym,msym+1:2*msym))))
% 	max(max(abs(aaa(1:msym,msym+1:2*msym))))
% 	max(max(abs(aaa(msym+1:2*msym,1:msym))))
% 	figure;spy(a);stop
if sym(1)~=0;  % symetries
alpha=sym(1);% alpha=sym(1)*exp(i*beta0*d);ATTENTION COMPRENDRE POURQUOI IL FAUT PAS BETA0  !!!

% si_elf=sparse([1:m-(mh+mb)].',[[1:msym].';[1:m-(mh+mb)-msym].'],[ones(msym,1);alpha*ones(msym,1);ones(m-(mh+mb)-2*msym,1)],m-(mh+mb),m-(mh+mb)-msym,m-(mh+mb));
% ssi_elf=sparse([[1:msym].';[1:m-(mh+mb)-msym].'],[1:m-(mh+mb)].',[.5*ones(msym,1);(.5/alpha)*ones(msym,1);ones(m-(mh+mb)-2*msym,1)],m-(mh+mb)-msym,m-(mh+mb),m-(mh+mb));

si_elf=sparse([1:m-(mh+mb)].',[[1:msym].';[1:m-(mh+mb)-msym].'],[ones(msym,1);alpha*ones(msym,1);ones(m-(mh+mb)-2*msym,1)],m-(mh+mb),m-(mh+mb)-msym,m-(mh+mb));
ssi_elf=sparse([[1:msym].';[1:m-(mh+mb)-msym].'],[1:m-(mh+mb)].',[.5*ones(msym,1);(.5/alpha)*ones(msym,1);ones(m-(mh+mb)-2*msym,1)],m-(mh+mb)-msym,m-(mh+mb),m-(mh+mb));


[l,u,p,q]=lu(ssi_elf*a(1:m-(mh+mb),1:m-(mh+mb))*si_elf);
%aa=a(m-(mh+mb)+1:m,m-(mh+mb)+1:m)-((a(m-(mh+mb)+1:m,1:m-(mh+mb))*(si_elf*q))/u)*(l\(p*(ssi_elf*a(1:m-(mh+mb),m-(mh+mb)+1:m))));aa=inv(full(aa));
aa=full(a(m-(mh+mb)+1:m,m-(mh+mb)+1:m))-retfastprod(((a(m-(mh+mb)+1:m,1:m-(mh+mb))*(si_elf*q))/u),(l\(p*(ssi_elf*a(1:m-(mh+mb),m-(mh+mb)+1:m)))));aa=inv(aa);

else
si_elf=[];ssi_elf=[];alpha=1;

[l,u,p,q]=lu(a(1:m-(mh+mb),1:m-(mh+mb)));
aa=full(a(m-(mh+mb)+1:m,m-(mh+mb)+1:m))-retfastprod(((a(m-(mh+mb)+1:m,1:m-(mh+mb))*q)/u),(l\(p*a(1:m-(mh+mb),m-(mh+mb)+1:m))));aa=inv(full(aa));

%aaa=a(m-(mh+mb)+1:m,m-(mh+mb)+1:m)-((a(m-(mh+mb)+1:m,1:m-(mh+mb))*q)/u)*(l\(p*a(1:m-(mh+mb),m-(mh+mb)+1:m)));aa=inv(full(aa));

end;

clear l u p q;

a_elf={a,aa,mh,mb,m,me,Hx,Hy,Bx,grad_grad,x,maille,indices,ep,si_elf,ssi_elf,alpha,struct('dim',1,'genre',2,'type',4,'sog',sog)};


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%                  ELEMENTS FINIS    2 D                 %%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [a_elf,m,me,maillage]=retelf_couche_2D(init,w,parm,k0,maillage);
if ~isfield(parm,'c');parm.c=0;end;if ~isfield(parm,'pol');parm.pol=0;end;
cao=init{end}.cao;d=init{end}.d;beta0=init{end}.beta0;sym=init{end}.sym;sog=init{end}.sog;if isempty(sym);sym=[0,0];end;
if isempty(cao);cao=[d/2,0,0];end;
if isempty(maillage);maillage=retmaillage(parm,w,d,sym,cao);end;  %  maillage
xcao=mod(cao(1),d(1));ycao=mod(cao(2),d(2));lcaox=real(cao(3));lcaoy=real(cao(4));caoix=imag(cao(3));caoiy=imag(cao(4));
% separation de la periode 2 D en  domaines pour le changement de coordonnees
if lcaox>0;x=mod([0,xcao+lcaox/2,xcao-lcaox/2],d(1));x=sort(retelimine([x,d(1)]));else x=[0,d(1)];end;nx=length(x)-1;
bornesx=[x(1:end-1);x(2:end)];centres=sum(bornesx,1)/2;
xc=xcao*ones(1,nx);
f=find(abs(xcao-centres)>(d(1)/2));xc(f)=xcao+d(1);
lx=lcaox*(  abs(centres-xc)<lcaox/2);

if lcaoy>0;y=mod([0,ycao+lcaoy/2,ycao-lcaoy/2],d(2));y=sort(retelimine([y,d(2)]));else y=[0,d(2)];end;ny=length(y)-1;
bornesy=[y(1:end-1);y(2:end)];centres=sum(bornesy,1)/2;
yc=ycao*ones(1,ny);
f=find(abs(ycao-centres)>(d(2)/2));yc(f)=ycao+d(2);
ly=lcaoy*(  abs(centres-yc)<lcaoy/2);
for ii=1:nx;for jj=1:ny;
domaines(ii+(jj-1)*nx)=struct('bornesx',bornesx(:,ii),'lx',lx(ii),'xc',xc(ii),'gx',1-1/(1+i*caoix),'gy',1-1/(1+i*caoiy),...
    'bornesy',bornesy(:,jj),'ly',ly(jj),'yc',yc(jj),'f',[],'Bx',[],'By',[],'bx',1,'by',1,...
'tf_tetra',[],'du_tf_tetra',[],'dv_tf_tetra',[],'dw_tf_tetra',[],'dudu_tf_tetra',[],'dvdv_tf_tetra',[],'dwdw_tf_tetra',[],'dvdw_tf_tetra',[],'dudw_tf_tetra',[]);
end;end;

%if isempty(maillage);maillage=retmaillage(parm,w,d,sym,cao);end;  %  maillage
ep=maillage.ep*k0;if parm.pol==2;ep=ep(:,[4,5,6,1,2,3]);end;% choix entre E (pol=0) et H (pol=2)
x=maillage.x;me=maillage.centre.me;m=maillage.centre.m;mh=maillage.haut.m;mb=maillage.bas.m;
indices=maillage.centre.indices;noeuds=maillage.centre.noeuds;arretes_local=maillage.centre.arretes_local;
aretes=maillage.centre.aretes;type=maillage.centre.type_arrete;
fac=[exp(i*beta0(1).*d(1)),exp(i*beta0(2).*d(2)),exp(i*beta0(1).*d(1)+i*beta0(2).*d(2))];
%%%%%%%%%%%%%%%%%%%%%%
%   ELEMENTS FINIS   %
%%%%%%%%%%%%%%%%%%%%%%

% volume des tetraedres  et fonctions w
grad3=ones(4,4,me);v=zeros(me,1);
for ii=1:4;for jj=1:3;grad3(ii,jj,:)=x(noeuds(:,ii),jj).';end;end;
for ii=1:me;grad3(:,:,ii)=inv(grad3(:,:,ii)).';v(ii)=1/(6*abs(det(grad3(1:3,1:3,ii))));end;


centres=(x(maillage.centre.noeuds(:,1),1:2)+x(maillage.centre.noeuds(:,2),1:2)+x(maillage.centre.noeuds(:,3),1:2)+x(maillage.centre.noeuds(:,4),1:2))/4;
for ii=1:length(domaines);% repartition en domaines ( au maximum 9)
domaines(ii).f=find((centres(:,1)>domaines(ii).bornesx(1))&(centres(:,1)<domaines(ii).bornesx(2))&(centres(:,2)>domaines(ii).bornesy(1))&(centres(:,2)<domaines(ii).bornesy(2)));
domaines(ii)=cala_init(domaines(ii),v(domaines(ii).f),maillage.centre.noeuds(domaines(ii).f,:),x);% calcul des matrices et des tf des tetraedres
end;
clear centres;
%  CONSTRUCTION DE LA MATRICE CREUSE DES ELEMENTS FINIS
aa=zeros(me,6,6);
for ii=1:6;
if (lcaox==0)&(lcaoy==0);jj0=ii;else jj0=1;end;% si pas cao matrice symetrique
for jj=jj0:6;
gi=zeros(me,4);gj=zeros(me,4);gii=zeros(me,4);gjj=zeros(me,4);
for k=1:4;% composantes x y z de grad(w)
gi(:,k)=grad3(   sub2ind([4,4,me],arretes_local(:,jj,1),k*ones(me,1),[1:me].') );      
gj(:,k)=grad3(   sub2ind([4,4,me],arretes_local(:,jj,2),k*ones(me,1),[1:me].') );      
gii(:,k)=grad3(  sub2ind([4,4,me],arretes_local(:,ii,1),k*ones(me,1),[1:me].') );      
gjj(:,k)=grad3(  sub2ind([4,4,me],arretes_local(:,ii,2),k*ones(me,1),[1:me].') );
end;  
aa(:,ii,jj)=v.*((gi(:,1).*gjj(:,1).*ep(indices,4)+gi(:,2).*gjj(:,2).*ep(indices,5)+gi(:,3).*gjj(:,3).*ep(indices,6)).*(.05*(1+ (arretes_local(:,jj,2)==arretes_local(:,ii,1)) ))... 
+(gj(:,1).*gii(:,1).*ep(indices,4)+gj(:,2).*gii(:,2).*ep(indices,5)+gj(:,3).*gii(:,3).*ep(indices,6)).*(.05*(1+ (arretes_local(:,jj,1)==arretes_local(:,ii,2)) ))... 
-(gi(:,1).*gii(:,1).*ep(indices,4)+gi(:,2).*gii(:,2).*ep(indices,5)+gi(:,3).*gii(:,3).*ep(indices,6)).*(.05*(1+ (arretes_local(:,jj,2)==arretes_local(:,ii,2)) ))... 
-(gj(:,1).*gjj(:,1).*ep(indices,4)+gj(:,2).*gjj(:,2).*ep(indices,5)+gj(:,3).*gjj(:,3).*ep(indices,6)).*(.05*(1+ (arretes_local(:,jj,1)==arretes_local(:,ii,1)) )));

for kk=1:length(domaines);f=domaines(kk).f;
aa(f,ii,jj)=aa(f,ii,jj)+cala(domaines(kk),gi(f,:),gj(f,:),gii(f,:),gjj(f,:),ep(indices(f),:),v(f));
end;
if (lcaox==0)&(lcaoy==0);if jj>ii;aa(:,jj,ii)=aa(:,ii,jj);end;end;% si pas cao matrice symetrique
end;
end;
aa=permute(aa,[2,3,1]);
% pseudo periodicite
for kk=1:3;
for ii=find(any(type==kk,2)).';
k=type(ii,:);f=find(k==kk);
aaa=ones(1,6);aaa(f)=fac(kk);
aa(:,:,ii)=aa(:,:,ii).*((1./aaa.')*aaa);
end;
end;

% construction de la matrice a
k1=zeros(6,6,me);k2=zeros(6,6,me);
for ii=1:6;prv=repmat(aretes(:,ii).',6,1);
k1(ii,:,:)=reshape(prv,1,6,me);
k2(:,ii,:)=reshape(prv,6,1,me);
end;
a=sparse(k1(:),k2(:),aa(:),m,m,36*me);

% construction des domaines du haut et du bas pour les conditions aux limites
domaines=rmfield(domaines,{'Bx','By','tf_tetra','du_tf_tetra','dv_tf_tetra','dw_tf_tetra','dudu_tf_tetra','dvdv_tf_tetra','dwdw_tf_tetra','dvdw_tf_tetra','dudw_tf_tetra'});
domaines_h=domaines;% domaines en haut
centres=(x(maillage.haut.noeuds(:,1),1:2)+x(maillage.haut.noeuds(:,2),1:2)+x(maillage.haut.noeuds(:,3),1:2))/3;
for ii=1:length(domaines); 
domaines_h(ii).f=find((centres(:,1)>domaines(ii).bornesx(1))&(centres(:,1)<domaines(ii).bornesx(2))&(centres(:,2)>domaines(ii).bornesy(1))&(centres(:,2)<domaines(ii).bornesy(2)));
end;
domaines_b=domaines;% domaines en bas
centres=(x(maillage.bas.noeuds(:,1),1:2)+x(maillage.bas.noeuds(:,2),1:2)+x(maillage.bas.noeuds(:,3),1:2))/3;
for ii=1:length(domaines); 
domaines_b(ii).f=find((centres(:,1)>domaines(ii).bornesx(1))&(centres(:,1)<domaines(ii).bornesx(2))&(centres(:,2)>domaines(ii).bornesy(1))&(centres(:,2)<domaines(ii).bornesy(2)));
end;
clear domaines

% construction des matrices EH qui servent a calculer les champs
EH=cell(1,6);
if real(parm.c)==1;  %   <------------------
aa=zeros(me,6,6);
for ii=1:6;
gi=zeros(me,4);gj=zeros(me,4);
for k=1:4;% composantes x y z de grad(w)
gi(:,k)=grad3(   sub2ind([4,4,me],arretes_local(:,ii,1),k*ones(me,1),[1:me].') );      
gj(:,k)=grad3(   sub2ind([4,4,me],arretes_local(:,ii,2),k*ones(me,1),[1:me].') );      
end;
aa(:,ii,1)=gj(:,1).*gi(:,4)-gi(:,1).*gj(:,4);
aa(:,ii,2)=gj(:,2).*gi(:,4)-gi(:,2).*gj(:,4);
aa(:,ii,3)=gj(:,3).*gi(:,4)-gi(:,3).*gj(:,4);
aa(:,ii,4)=gj(:,3).*gi(:,2)-gi(:,3).*gj(:,2);
aa(:,ii,5)=gj(:,1).*gi(:,3)-gi(:,1).*gj(:,3);
aa(:,ii,6)=gj(:,2).*gi(:,1)-gi(:,2).*gj(:,1);
for kk=1:3;f=find(type(:,ii)==kk);aa(f,ii,:)=aa(f,ii,:)*fac(kk);end;% pseudo periodicite
end;
aa=reshape(permute(aa,[2,1,3]),6*me,6);
k1=repmat(1:me,6,1);k2=aretes.';
for ii=1:6;EH{ii}=sparse(k1(:),k2(:),aa(:,ii),me,m,6*me);end;
if ~isreal(parm.c);EH=retio(EH,1);end;
end;  %   <------------------
clear aa gi gii gj gjj k1 k2 prv v grad3;

%%%%%%%%w%%%%%%%%%%%%%%%%%%%
%  ELIMINATION  :aa*E=H   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%


% % 
% % 
%             %%%%%%%%		test de symetrie
%              figure;spy(a);
% 
%             if abs(sym(1))+abs(sym(2))==1;
%             mprv=round(maillage.nb_sym.nsym/2);
%             aaa=a(1:2*mprv,1:2*mprv);
%             max(max(abs(aaa(1:mprv,1:mprv)-aaa(mprv+1:2*mprv,mprv+1:2*mprv))))
%             max(max(abs(aaa(1:mprv,mprv+1:2*mprv)))),max(max(abs(aaa(mprv+1:2*mprv,1:mprv))))
%             end;       
% 
%             if abs(sym(1))+abs(sym(2))==2;
%             mprv=round(maillage.nb_sym.nsym/4);
%             aaa=a(1:4*mprv,1:4*mprv);
%             max(max(abs(aaa(1:mprv,1:mprv)-aaa(mprv+1:2*mprv,mprv+1:2*mprv))))
%             max(max(abs(aaa(1:mprv,1:mprv)-aaa(2*mprv+1:3*mprv,2*mprv+1:3*mprv))))
%             max(max(abs(aaa(1:mprv,1:mprv)-aaa(3*mprv+1:4*mprv,3*mprv+1:4*mprv))))
%             max(max(abs(aaa(1:mprv,mprv+1:4*mprv)))),max(max(abs(aaa(mprv+1:4*mprv,1:mprv))))
%             max(max(abs(aaa(1:mprv,mprv+1:4*mprv)))),max(max(abs(aaa(mprv+1:4*mprv,1:mprv))))
%             max(max(abs(aaa(mprv+1:2*mprv,2*mprv+1:4*mprv)))),max(max(abs(aaa(2*mprv+1:4*mprv,mprv+1:2*mprv))))
%             max(max(abs(aaa(2*mprv+1:3*mprv,3*mprv+1:4*mprv)))),max(max(abs(aaa(3*mprv+1:4*mprv,2*mprv+1:3*mprv))))
%             end;       
%             %stop
%             
clear indices noeuds arretes_local aretes type x;
if parm.pol==0;fac_sym=1;else;fac_sym=-1;end;
switch(abs(sym(1))+i*abs(sym(2)));% <----symetries
case 1;   % symetrie en x  
msym=maillage.nb_sym.nsym/2;
alphay=1;
%if parm.pol==0;alphax=-sym(1)*exp(i*beta0(1)*d(1));else;alphax=sym(1)*exp(i*beta0(1)*d(1));end;
if parm.pol==0;alphax=-sym(1);else;alphax=sym(1);end;% ATTENTION COMPRENDRE POURQUOI IL FAUT PAS BETA0  !!!

si_elf=sparse([1:m-mb-mh].',[1:msym,1:m-mb-mh-msym].',[(1/alphax)*ones(msym,1);ones(m-mb-mh-msym,1)],m-mb-mh,m-mb-mh-msym,m-mb-mh);
ssi_elf=sparse([1:m-mb-mh-msym].',[msym+1:m-mb-mh].',ones(m-mb-mh-msym,1),m-mb-mh-msym,m-mb-mh,m-mb-mh-msym);
msymh=maillage.nb_sym.nsymh/2;nnsymh=maillage.nb_sym.nnsymh;msymb=maillage.nb_sym.nsymb/2;nnsymb=maillage.nb_sym.nnsymb;
v1=(1:2*msymh+nnsymh+2*msymb+nnsymb).';
v2=[[1:msymh].';[1:msymh+nnsymh].';msymh+nnsymh+[[1:msymb].';[1:msymb+nnsymb].']];
uu=[(1/alphax)*ones(msymh,1);ones(msymh+nnsymh,1);(1/alphax)*ones(msymb,1);ones(msymb+nnsymb,1)];
si_bord=sparse(v1,v2,uu,2*msymh+nnsymh+2*msymb+nnsymb,msymh+nnsymh+msymb+nnsymb,2*msymh+nnsymh+2*msymb+nnsymb);
v1=[1:msymh+nnsymh+msymb+nnsymb].';
v2=[[msymh+1:2*msymh+nnsymh].';2*msymh+nnsymh+msymb+[1:msymb+nnsymb].'];
uu=ones(msymh+nnsymh+msymb+nnsymb,1);
ssi_bord=sparse(v1,v2,uu,msymh+nnsymh+msymb+nnsymb,2*msymh+nnsymh+2*msymb+nnsymb,msymh+nnsymh+msymb+nnsymb);

case i;   % symetrie en y 
msym=maillage.nb_sym.nsym/2;
alphax=1;
%if parm.pol==0;alphay=sym(2)*exp(i*beta0(2)*d(2));else;alphay=-sym(2)*exp(i*beta0(2)*d(2));end;% ATTENTION COMPRENDRE POURQUOI IL FAUT PAS BETA0  !!!
if parm.pol==0;alphay=sym(2);else;alphay=-sym(2);end;% ATTENTION COMPRENDRE POURQUOI IL FAUT PAS BETA0  !!!
% 
si_elf=sparse([1:m-mb-mh].',[1:msym,1:m-mb-mh-msym].',[(1/alphay)*ones(msym,1);ones(m-mb-mh-msym,1)],m-mb-mh,m-mb-mh-msym,m-mb-mh);
ssi_elf=sparse([1:m-mb-mh-msym].',[msym+1:m-mb-mh].',ones(m-mb-mh-msym,1),m-mb-mh-msym,m-mb-mh,m-mb-mh-msym);

msymh=maillage.nb_sym.nsymh/2;nnsymh=maillage.nb_sym.nnsymh;msymb=maillage.nb_sym.nsymb/2;nnsymb=maillage.nb_sym.nnsymb;
v1=(1:2*msymh+nnsymh+2*msymb+nnsymb).';
v2=[[1:msymh].';[1:msymh+nnsymh].';msymh+nnsymh+[[1:msymb].';[1:msymb+nnsymb].']];
uu=[(1/alphay)*ones(msymh,1);ones(msymh+nnsymh,1);(1/alphay)*ones(msymb,1);ones(msymb+nnsymb,1)];
si_bord=sparse(v1,v2,uu,2*msymh+nnsymh+2*msymb+nnsymb,msymh+nnsymh+msymb+nnsymb,2*msymh+nnsymh+2*msymb+nnsymb);
v1=[1:msymh+nnsymh+msymb+nnsymb].';
v2=[[msymh+1:2*msymh+nnsymh].';2*msymh+nnsymh+msymb+[1:msymb+nnsymb].'];
uu=ones(msymh+nnsymh+msymb+nnsymb,1);
ssi_bord=sparse(v1,v2,uu,msymh+nnsymh+msymb+nnsymb,2*msymh+nnsymh+2*msymb+nnsymb,msymh+nnsymh+msymb+nnsymb);

case 1+i; % symetrie en x et  y  
msym=maillage.nb_sym.nsym/4;

%if parm.pol==0;alphay=sym(2)*exp(i*beta0(2)*d(2));alphax=-sym(1)*exp(i*beta0(1)*d(1));else;alphay=-sym(2)*exp(i*beta0(2)*d(2));alphax=sym(1)*exp(i*beta0(1)*d(1));end;
if parm.pol==0;alphay=sym(2);alphax=-sym(1);else;alphay=-sym(2);alphax=sym(1);end;% ATTENTION COMPRENDRE POURQUOI IL FAUT PAS BETA0  !!!


si_elf=sparse([1:m-mb-mh].',[1:msym,1:msym,1:msym,1:m-mb-mh-3*msym].',[1/(alphax*alphay)*ones(msym,1);(1/alphay)*ones(msym,1);(1/alphax)*ones(msym,1);ones(m-mb-mh-3*msym,1)],m-mb-mh,m-mb-mh-3*msym,m-mb-mh);
ssi_elf=sparse([1:m-mb-mh-3*msym].',[3*msym+1:m-mb-mh].',ones(m-mb-mh-3*msym,1),m-mb-mh-3*msym,m-mb-mh,m-mb-mh-3*msym);

msymh=maillage.nb_sym.nsymh/4;nnsymh=maillage.nb_sym.nnsymh;msymb=maillage.nb_sym.nsymb/4;nnsymb=maillage.nb_sym.nnsymb;
v1=(1:4*msymh+nnsymh+4*msymb+nnsymb).';
v2=[[1:msymh].';[1:msymh].';[1:msymh].';[1:msymh+nnsymh].';msymh+nnsymh+[[1:msymb].';[1:msymb].';[1:msymb].';[1:msymb+nnsymb].']];
%uu=[alphax*ones(msymh,1);alphay*ones(msymh,1);alphax*alphay*ones(msymh,1);ones(msymh+nnsymh,1);...
%alphax*ones(msymb,1);alphay*ones(msymb,1);alphax*alphay*ones(msymb,1);ones(msymb+nnsymb,1)];

uu=[1/(alphax*alphay)*ones(msymh,1);(1/alphay)*ones(msymh,1);(1/alphax)*ones(msymh,1);ones(msymh+nnsymh,1);...
1/(alphax*alphay)*ones(msymb,1);(1/alphay)*ones(msymb,1);(1/alphax)*ones(msymb,1);ones(msymb+nnsymb,1)];
si_bord=sparse(v1,v2,uu,4*msymh+nnsymh+4*msymb+nnsymb,msymh+nnsymh+msymb+nnsymb,4*msymh+nnsymh+4*msymb+nnsymb);
v1=[1:msymh+nnsymh+msymb+nnsymb].';
v2=[[3*msymh+1:4*msymh+nnsymh].';4*msymh+nnsymh+3*msymb+[1:msymb+nnsymb].'];
uu=ones(msymh+nnsymh+msymb+nnsymb,1);
ssi_bord=sparse(v1,v2,uu,msymh+nnsymh+msymb+nnsymb,4*msymh+nnsymh+4*msymb+nnsymb,msymh+nnsymh+msymb+nnsymb);

case 0;   % pas de symetrie  
alphax=1;alphay=1;si_elf=[];ssi_elf=[];si_bord=[];ssi_bord=[];
end;%  <-------symetries
clear v1 v2;

if isempty(si_elf);
[l,u,p,q]=lu(a(1:m-mb-mh,1:m-mb-mh));
aa=a(m-mb-mh+1:m,1:m-mb-mh)*q;clear q;aa=(u.'\aa.').';clear u;
aaa=p*a(1:m-mb-mh,m-mb-mh+1:m);clear p;aaa=(l\aaa);clear l;
aa=a(m-mb-mh+1:m,m-mb-mh+1:m)-aa*aaa;aa=inv(full(aa));
else;
[l,u,p,q]=lu(ssi_elf*a(1:m-mb-mh,1:m-mb-mh)*si_elf);
aa=((ssi_bord*(a(m-mb-mh+1:m,1:m-mb-mh))*si_elf)*q);clear q;aa=(u.'\aa.').';clear u;
aaa=(p*(ssi_elf*(a(1:m-mb-mh,m-mb-mh+1:m)*si_bord)));clear p;aaa=(l\aaa);clear l;
saa=size(aa);

if abs(sym(1))+i*abs(sym(2))==1+i; % 2 symetries  pour gain de temps on decompose le produit aa*aa en sparse et full
aa=full(ssi_bord*a(m-mb-mh+1:m,m-mb-mh+1:m)*si_bord-aa(:,1:saa(2)-saa(1)-1)*aaa(1:saa(2)-saa(1)-1,:))...
-full(aa(:,saa(2)-saa(1):saa(2)))*full(aaa(saa(2)-saa(1):saa(2),:));
else;
aa=ssi_bord*a(m-mb-mh+1:m,m-mb-mh+1:m)*si_bord-aa*aaa;
end;    
clear aaa;aa=inv(full(aa));
end;

% %  construction des matrices qui permettent la liaison a RCWA
% 	Ah1,Bh1,Ch1 font passer des aretes aux triangles dans la pseudo periodicite on divise par fac
% 	Ah2,Bh2,Ch2 font passer des  triangles aux aretes dans la pseudo periodicite on multiplie par fac
%    idem en bas pour Ab1,Bb1,Cb1 ,Ab2,Bb2,Cb2 
me=maillage.centre.me;m=maillage.centre.m;mh=maillage.haut.m;mb=maillage.bas.m;
[Ah1,Bh1,Ch1]=calABC(maillage.haut,maillage.x,m-maillage.haut.m-maillage.bas.m,1./fac); % en haut
[Ab1,Bb1,Cb1]=calABC(maillage.bas,maillage.x,m-maillage.bas.m,1./fac);                 % en bas             
if all(beta0==0);[Ah2,Bh2,Ch2,Ab2,Bb2,Cb2]=deal(Ah1,Bh1,Ch1,Ab1,Bb1,Cb1);
else;
[Ah2,Bh2,Ch2]=calABC(maillage.haut,maillage.x,m-maillage.haut.m-maillage.bas.m,fac); % en haut
[Ab2,Bb2,Cb2]=calABC(maillage.bas,maillage.x,m-maillage.bas.m,fac);                 % en bas
end;
if ~isempty(si_elf);
Ah2=Ah2*si_bord(1:mh,1:msymh+nnsymh);
Bh2=Bh2*si_bord(1:mh,1:msymh+nnsymh);
Ch2=Ch2*si_bord(1:mh,1:msymh+nnsymh);
Ah1=Ah1*ssi_bord(1:msymh+nnsymh,1:mh).';
Bh1=Bh1*ssi_bord(1:msymh+nnsymh,1:mh).';
Ch1=Ch1*ssi_bord(1:msymh+nnsymh,1:mh).';

Ab2=Ab2*si_bord(mh+1:end,msymh+nnsymh+1:end);
Bb2=Bb2*si_bord(mh+1:end,msymh+nnsymh+1:end);
Cb2=Cb2*si_bord(mh+1:end,msymh+nnsymh+1:end);
Ab1=Ab1*ssi_bord(msymh+nnsymh+1:end,mh+1:end).';
Bb1=Bb1*ssi_bord(msymh+nnsymh+1:end,mh+1:end).';
Cb1=Cb1*ssi_bord(msymh+nnsymh+1:end,mh+1:end).';

end;

a_elf={a,aa,maillage,Ah1.',Ah2,Ab1.',Ab2,Bh1.',Bh2,Bb1.',Bb2,Ch1.',Ch2,Cb1.',Cb2,si_elf,ssi_elf,si_bord,ssi_bord,alphax,alphay,EH,ep,domaines_h,domaines_b,struct('dim',1,'genre',2,'type',4,'sog',sog,'pol',parm.pol)};

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [A,B,C]=calABC(maillage,x,offset,fac);
mme=maillage.me;
aarretes=maillage.aretes-offset;
%  fonctions w en 2D
grad2=ones(3,3,mme);AA=zeros(mme,3);BB=zeros(mme,3);CC=zeros(mme,3);
for ii=1:3;for jj=1:2;grad2(ii,jj,:)=x(maillage.noeuds(:,ii),jj).';end;end;
for ii=1:mme;grad2(:,:,ii)=inv(grad2(:,:,ii)).';end;
for ii=1:3;% <------ 
gi=zeros(mme,3);gj=zeros(mme,3);
for k=1:3;% composantes x y z de grad(w)
gi(:,k)=grad2(   sub2ind([3,3,mme],maillage.arretes_local(:,ii,1),k*ones(mme,1),[1:mme].') );      
gj(:,k)=grad2(   sub2ind([3,3,mme],maillage.arretes_local(:,ii,2),k*ones(mme,1),[1:mme].') );      
end;  
AA(:,ii)=gj(:,1).*gi(:,2)-gi(:,1).*gj(:,2);
BB(:,ii)=gi(:,2).*gj(:,3)-gj(:,2).*gi(:,3);
CC(:,ii)=gj(:,1).*gi(:,3)-gi(:,1).*gj(:,3);
% pseudo periodicite (pas de type 3 )
for kk=1:2;f=find(maillage.type_arrete(:,ii)==kk);AA(f,ii)=AA(f,ii)*fac(kk);BB(f,ii)=BB(f,ii)*fac(kk);CC(f,ii)=CC(f,ii)*fac(kk);end;
end;     %  <------
A=sparse(repmat([1:mme]',3,1),aarretes(:),AA(:),mme,maillage.m,3*mme);
B=sparse(repmat([1:mme]',3,1),aarretes(:),BB(:),mme,maillage.m,3*mme);
C=sparse(repmat([1:mme]',3,1),aarretes(:),CC(:),mme,maillage.m,3*mme);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function domaines=cala_init(domaines,v,noeuds,x);
sz=size(noeuds);xx=reshape(x(noeuds,1),sz);yy=reshape(x(noeuds,2),sz);zz=reshape(x(noeuds,3),sz);delta=1.e-4;
xc=sum(xx,2)/4;yc=sum(yy,2)/4;zc=sum(zz,2)/4;xx=xx-repmat(xc,1,4);yy=yy-repmat(yc,1,4);zz=zz-repmat(zc,1,4);% recentrement
switch((domaines.lx>0)+i*(domaines.ly>0));
case 0;return;% pas de changement de coordonnees
case 1;%  changement de coordonnees en x seulement
% calcul de Bx
alpx=2i*pi*domaines.xc/domaines.lx;
Bx=[exp(-2*alpx)*domaines.gx/16,-.25*exp(-alpx),.5-domaines.gx/8,-.25*exp(alpx),exp(2*alpx)*domaines.gx/16];
%    x=linspace(domaines.bornesx(1),domaines.bornesx(2),101);phi=pi*(x-domaines.xc)/domaines.lx;fx=(1-cos(phi).^2).*(1-domaines.gx*(cos(phi)).^2);    
%     prv=zeros(size(x));for ii=-2:2;prv=prv+Bx(ii+3)*exp(-2i*pi*ii*x/domaines.lx);end;retcompare(prv,(fx))
%     figure;plot(x,real(fx),x,imag(fx));
domaines.bx=Bx;Bx=conv(Bx,Bx);domaines.Bx=Bx;
%       prv=zeros(size(x));for ii=-4:4;prv=prv+Bx(ii+5)*exp(-2i*pi*ii*x/domaines.lx);end;retcompare(prv,fx.^2)
ux=[-4:4].'/domaines.lx;uy=zeros(9,1);uz=zeros(9,1);

[domaines.tf_tetra,prv]=rettfsimplex(xx.',yy.',zz.',ux,uy,uz,v);
domaines.du_tf_tetra=prv(:,:,1,1);domaines.dv_tf_tetra=prv(:,:,1,2);domaines.dw_tf_tetra=prv(:,:,1,3);
[prv,prv_p]=rettfsimplex(xx.',yy.',zz.',ux+delta,uy,uz,v);[prv,prv_m]=rettfsimplex(xx.',yy.',zz.',ux-delta,uy,uz,v);
domaines.dudu_tf_tetra=(prv_p(:,:,1,1)-prv_m(:,:,1,1))/(2*delta);
[prv,prv_p]=rettfsimplex(xx.',yy.',zz.',ux,uy+delta,uz,v);[prv,prv_m]=rettfsimplex(xx.',yy.',zz.',ux,uy-delta,uz,v);
domaines.dvdv_tf_tetra=(prv_p(:,:,1,2)-prv_m(:,:,1,2))/(2*delta);
[prv,prv_p]=rettfsimplex(xx.',yy.',zz.',ux,uy,uz+delta,v);[prv,prv_m]=rettfsimplex(xx.',yy.',zz.',ux,uy,uz-delta,v);
domaines.dwdw_tf_tetra=(prv_p(:,:,1,3)-prv_m(:,:,1,3))/(2*delta);
domaines.dvdw_tf_tetra=(prv_p(:,:,1,2)-prv_m(:,:,1,2))/(2*delta);
m=9;


case i;%  changement de coordonnees en y seulement
% calcul de By
alpy=2i*pi*domaines.yc/domaines.ly;
By=[exp(-2*alpy)*domaines.gy/16,-.25*exp(-alpy),.5-domaines.gy/8,-.25*exp(alpy),exp(2*alpy)*domaines.gy/16];
%        y=linspace(domaines.bornesy(1),domaines.bornesy(2),101);phi=pi*(y-domaines.yc)/domaines.ly;fy=(1-cos(phi).^2).*(1-domaines.gy*(cos(phi)).^2);    
%       prv=zeros(size(y));for ii=-2:2;prv=prv+By(ii+3)*exp(-2i*pi*ii*y/domaines.ly);end;retcompare(prv,(fy))
%       figure;plot(y,real(fy),y,imag(fy));
domaines.by=By;By=conv(By,By);domaines.By=By;
uy=[-4:4].'/domaines.ly;ux=zeros(9,1);uz=zeros(9,1);


[domaines.tf_tetra,prv]=rettfsimplex(xx.',yy.',zz.',ux,uy,uz,v);
domaines.du_tf_tetra=prv(:,:,1,1);domaines.dv_tf_tetra=prv(:,:,1,2);domaines.dw_tf_tetra=prv(:,:,1,3);
[prv,prv_p]=rettfsimplex(xx.',yy.',zz.',ux+delta,uy,uz,v);[prv,prv_m]=rettfsimplex(xx.',yy.',zz.',ux-delta,uy,uz,v);
domaines.dudu_tf_tetra=(prv_p(:,:,1,1)-prv_m(:,:,1,1))/(2*delta);
[prv,prv_p]=rettfsimplex(xx.',yy.',zz.',ux,uy+delta,uz,v);[prv,prv_m]=rettfsimplex(xx.',yy.',zz.',ux,uy-delta,uz,v);
domaines.dvdv_tf_tetra=(prv_p(:,:,1,2)-prv_m(:,:,1,2))/(2*delta);
[prv,prv_p]=rettfsimplex(xx.',yy.',zz.',ux,uy,uz+delta,v);[prv,prv_m]=rettfsimplex(xx.',yy.',zz.',ux,uy,uz-delta,v);
domaines.dwdw_tf_tetra=(prv_p(:,:,1,3)-prv_m(:,:,1,3))/(2*delta);
domaines.dudw_tf_tetra=(prv_p(:,:,1,1)-prv_m(:,:,1,1))/(2*delta);
m=9;
%        prv=zeros(size(y));for ii=-4:4;prv=prv+By(ii+5)*exp(-2i*pi*ii*y/domaines.ly);end;retcompare(prv,fy.^2)
case 1+i;%  changement de coordonnees en x  et y
% calcul de Bx  fx^2
alpx=2i*pi*domaines.xc/domaines.lx;
Bx=[exp(-2*alpx)*domaines.gx/16,-.25*exp(-alpx),.5-domaines.gx/8,-.25*exp(alpx),exp(2*alpx)*domaines.gx/16];
domaines.bx=Bx;Bx=conv(Bx,Bx);domaines.Bx=Bx;
alpy=2i*pi*domaines.yc/domaines.ly;
By=[exp(-2*alpy)*domaines.gy/16,-.25*exp(-alpy),.5-domaines.gy/8,-.25*exp(alpy),exp(2*alpy)*domaines.gy/16];
domaines.by=By;By=conv(By,By);domaines.By=By;

ux=[-4:4].'/domaines.lx;uy=[-4:4].'/domaines.ly;[ux,uy]=retmeshgrid(ux,uy);uz=zeros(size(ux));
[domaines.tf_tetra,prv]=rettfsimplex(xx.',yy.',zz.',ux,uy,uz,v);
domaines.du_tf_tetra=prv(:,:,:,1);domaines.dv_tf_tetra=prv(:,:,:,2);domaines.dw_tf_tetra=prv(:,:,:,3);
[prv,prv_p]=rettfsimplex(xx.',yy.',zz.',ux+delta,uy,uz,v);[prv,prv_m]=rettfsimplex(xx.',yy.',zz.',ux-delta,uy,uz,v);
domaines.dudu_tf_tetra=(prv_p(:,:,:,1)-prv_m(:,:,:,1))/(2*delta);
[prv,prv_p]=rettfsimplex(xx.',yy.',zz.',ux,uy+delta,uz,v);[prv,prv_m]=rettfsimplex(xx.',yy.',zz.',ux,uy-delta,uz,v);
domaines.dvdv_tf_tetra=(prv_p(:,:,:,2)-prv_m(:,:,:,2))/(2*delta);
[prv,prv_p]=rettfsimplex(xx.',yy.',zz.',ux,uy,uz+delta,v);[prv,prv_m]=rettfsimplex(xx.',yy.',zz.',ux,uy,uz-delta,v);
domaines.dwdw_tf_tetra=(prv_p(:,:,:,3)-prv_m(:,:,:,3))/(2*delta);
domaines.dudw_tf_tetra=(prv_p(:,:,:,1)-prv_m(:,:,:,1))/(2*delta);
domaines.dvdw_tf_tetra=(prv_p(:,:,:,2)-prv_m(:,:,:,2))/(2*delta);
m=81;

end;

% correction des derivees a cause du recentrement
exp_c=exp(-2i*pi*(xc(:)*ux(:).'+yc(:)*uy(:).'+zc(:)*uz(:).'));xc=repmat(xc,1,m);yc=repmat(yc,1,m);zc=repmat(zc,1,m);
exp_c=exp_c(:);xc=xc(:);yc=yc(:);zc=zc(:);
domaines.tf_tetra(:)=exp_c.*domaines.tf_tetra(:);
domaines.du_tf_tetra(:)=exp_c.*domaines.du_tf_tetra(:);
domaines.dv_tf_tetra(:)=exp_c.*domaines.dv_tf_tetra(:); 
domaines.dw_tf_tetra(:)=exp_c.*domaines.dw_tf_tetra(:); 
domaines.dudu_tf_tetra(:)=exp_c.*domaines.dudu_tf_tetra(:); 
domaines.dvdv_tf_tetra(:)=exp_c.*domaines.dvdv_tf_tetra(:); 
domaines.dwdw_tf_tetra(:)=exp_c.*domaines.dwdw_tf_tetra(:); 
if ~isempty(domaines.dudw_tf_tetra);domaines.dudw_tf_tetra(:)=exp_c.*domaines.dudw_tf_tetra(:);end; 
if ~isempty(domaines.dvdw_tf_tetra);domaines.dvdw_tf_tetra(:)=exp_c.*domaines.dvdw_tf_tetra(:);end;  
domaines.dudu_tf_tetra(:)=-4*pi^2*(xc.^2).*domaines.tf_tetra(:)-4i*pi*xc.*domaines.du_tf_tetra(:)+domaines.dudu_tf_tetra(:);
domaines.dvdv_tf_tetra(:)=-4*pi^2*(yc.^2).*domaines.tf_tetra(:)-4i*pi*yc.*domaines.dv_tf_tetra(:)+domaines.dvdv_tf_tetra(:);
domaines.dwdw_tf_tetra(:)=-4*pi^2*(zc.^2).*domaines.tf_tetra(:)-4i*pi*zc.*domaines.dw_tf_tetra(:)+domaines.dwdw_tf_tetra(:);
if ~isempty(domaines.dudw_tf_tetra);domaines.dudw_tf_tetra(:)=-4*pi^2*xc.*zc.*domaines.tf_tetra(:)-2i*pi*zc.*domaines.du_tf_tetra(:)-2i*pi*xc.*domaines.dw_tf_tetra(:)+domaines.dudw_tf_tetra(:);end;
if ~isempty(domaines.dvdw_tf_tetra);domaines.dvdw_tf_tetra(:)=-4*pi^2*yc.*zc.*domaines.tf_tetra(:)-2i*pi*zc.*domaines.dv_tf_tetra(:)-2i*pi*yc.*domaines.dw_tf_tetra(:)+domaines.dvdw_tf_tetra(:);end;
domaines.du_tf_tetra(:)=-2i*pi*xc.*domaines.tf_tetra(:)+domaines.du_tf_tetra(:);
domaines.dv_tf_tetra(:)=-2i*pi*yc.*domaines.tf_tetra(:)+domaines.dv_tf_tetra(:);
domaines.dw_tf_tetra(:)=-2i*pi*zc.*domaines.tf_tetra(:)+domaines.dw_tf_tetra(:);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function a=cala(domaines,gi,gj,gii,gjj,ep,v);

hx=gj(:,3).*gi(:,2)-gi(:,3).*gj(:,2);hy=gj(:,1).*gi(:,3)-gi(:,1).*gj(:,3);hz=gj(:,2).*gi(:,1)-gi(:,2).*gj(:,1);
hhx=gjj(:,3).*gii(:,2)-gii(:,3).*gjj(:,2);hhy=gjj(:,1).*gii(:,3)-gii(:,1).*gjj(:,3);hhz=gjj(:,2).*gii(:,1)-gii(:,2).*gjj(:,1);

switch((domaines.lx>0)+i*(domaines.ly>0));
    
case 0;%    pas de changement de coordonnees
hx=hx./ep(:,1);hy=hy./ep(:,2);hz=hz./ep(:,3);
a=4*v.*(hx.*hhx+hy.*hhy+hz.*hhz);

case 1;%    changement de coordonnees en x seulement
eex=gjj(:,1).*gii(:,4)-gii(:,1).*gjj(:,4);eey=gjj(:,2).*gii(:,4)-gii(:,2).*gjj(:,4);eez=gjj(:,3).*gii(:,4)-gii(:,3).*gjj(:,4);
ex=gj(:,1).*gi(:,4)-gi(:,1).*gj(:,4);
a=zeros(size(v));

BBx=-domaines.Bx;BBx(5)=BBx(5)+1;

for ii=-4:4;
a=a+BBx(ii+5)*ep(:,4).*(       ex.*eex.*domaines.tf_tetra(:,ii+5)...
+((hz.*eex+hhz.*ex).*domaines.dv_tf_tetra(:,ii+5)-(hy.*eex+hhy.*ex).*domaines.dw_tf_tetra(:,ii+5))/(2i*pi)...
+ ((hz.*hhy+hhz.*hy).*domaines.dvdw_tf_tetra(:,ii+5)-hz.*hhz.*domaines.dvdv_tf_tetra(:,ii+5)-hy.*hhy.*domaines.dwdw_tf_tetra(:,ii+5))/(4*(pi^2))); 
end;


hx=hx./ep(:,1);hy=hy./ep(:,2);hz=hz./ep(:,3);
a=a+4*v.*hx.*hhx;

for ii=-4:4;
a=a+domaines.Bx(ii+5)*(hy.*((4*hhy+(2i*pi*ii/domaines.lx)*eez).*domaines.tf_tetra(:,ii+5)...
+(ii/domaines.lx)*(hhy.*domaines.du_tf_tetra(:,ii+5)-hhx.*domaines.dv_tf_tetra(:,ii+5)))...
+hz.*((4*hhz-(2i*pi*ii/domaines.lx)*eey).*domaines.tf_tetra(:,ii+5)...
+(ii/domaines.lx).*(hhz.*domaines.du_tf_tetra(:,ii+5)-hhx.*domaines.dw_tf_tetra(:,ii+5)))) ;

end;

case i;%    changement de coordonnees en y seulement
eex=gjj(:,1).*gii(:,4)-gii(:,1).*gjj(:,4);eey=gjj(:,2).*gii(:,4)-gii(:,2).*gjj(:,4);eez=gjj(:,3).*gii(:,4)-gii(:,3).*gjj(:,4);
ey=gj(:,2).*gi(:,4)-gi(:,2).*gj(:,4);
a=zeros(size(v));

BBy=-domaines.By;BBy(5)=BBy(5)+1;
for ii=-4:4;
a=a+BBy(ii+5)*ep(:,5).*(       ey.*eey.*domaines.tf_tetra(:,ii+5)...
-((hz.*eey+hhz.*ey).*domaines.du_tf_tetra(:,ii+5)-(hx.*eey+hhx.*ey).*domaines.dw_tf_tetra(:,ii+5))/(2i*pi)...
+((hz.*hhx+hhz.*hx).*domaines.dudw_tf_tetra(:,ii+5)-hz.*hhz.*domaines.dudu_tf_tetra(:,ii+5)-hx.*hhx.*domaines.dwdw_tf_tetra(:,ii+5))/(4*pi^2));   
end;

hx=hx./ep(:,1);hy=hy./ep(:,2);hz=hz./ep(:,3);
a=a+4*v.*hy.*hhy;

for ii=-4:4;
a=a+domaines.By(ii+5)*(hx.*((4*hhx-(2i*pi*ii/domaines.ly)*eez).*domaines.tf_tetra(:,ii+5)...
+(ii/domaines.ly)*(-hhy.*domaines.du_tf_tetra(:,ii+5)+hhx.*domaines.dv_tf_tetra(:,ii+5)))...
+hz.*((4*hhz+(2i*pi*ii/domaines.ly)*eex).*domaines.tf_tetra(:,ii+5)...
+(ii/domaines.ly).*(hhz.*domaines.dv_tf_tetra(:,ii+5)-hhy.*domaines.dw_tf_tetra(:,ii+5))));
end;

case 1+i;%  changement de coordonnees en x  et y
ex=gj(:,1).*gi(:,4)-gi(:,1).*gj(:,4);ey=gj(:,2).*gi(:,4)-gi(:,2).*gj(:,4);ez=gj(:,3).*gi(:,4)-gi(:,3).*gj(:,4);
eex=gjj(:,1).*gii(:,4)-gii(:,1).*gjj(:,4);eey=gjj(:,2).*gii(:,4)-gii(:,2).*gjj(:,4);eez=gjj(:,3).*gii(:,4)-gii(:,3).*gjj(:,4);
a=zeros(size(v));

BBx=-domaines.Bx;BBx(5)=BBx(5)+1;
for ii=-4:4;
a=a+BBx(ii+5)*ep(:,4).*(       ex.*eex.*domaines.tf_tetra(:,ii+5,5)...
+((hz.*eex+hhz.*ex).*domaines.dv_tf_tetra(:,ii+5,5)-(hy.*eex+hhy.*ex).*domaines.dw_tf_tetra(:,ii+5,5))/(2i*pi)...
+ ((hz.*hhy+hhz.*hy).*domaines.dvdw_tf_tetra(:,ii+5,5)-hz.*hhz.*domaines.dvdv_tf_tetra(:,ii+5,5)-hy.*hhy.*domaines.dwdw_tf_tetra(:,ii+5,5))/(4*pi^2)); 
end;

BBy=-domaines.By;BBy(5)=BBy(5)+1;
for ii=-4:4;
a=a+BBy(ii+5)*ep(:,5).*(     ey.*eey.*domaines.tf_tetra(:,5,ii+5)...
-((hz.*eey+hhz.*ey).*domaines.du_tf_tetra(:,5,ii+5)-(hx.*eey+hhx.*ey).*domaines.dw_tf_tetra(:,5,ii+5))/(2i*pi)...
+ ((hz.*hhx+hhz.*hx).*domaines.dudw_tf_tetra(:,5,ii+5)-hz.*hhz.*domaines.dudu_tf_tetra(:,5,ii+5)-hx.*hhx.*domaines.dwdw_tf_tetra(:,5,ii+5))/(4*pi^2));    
end;

hx=hx./ep(:,1);hy=hy./ep(:,2);hz=hz./ep(:,3);
for ii=-4:4;
a=a+domaines.Bx(ii+5)*hy.*((4*hhy+(2i*pi*ii/domaines.lx)*eez).*domaines.tf_tetra(:,ii+5,5)...
+(ii/domaines.lx)*(hhy.*domaines.du_tf_tetra(:,ii+5,5)-hhx.*domaines.dv_tf_tetra(:,ii+5,5)))...
+domaines.By(ii+5)*hx.*((4*hhx-(2i*pi*ii/domaines.ly)*eez).*domaines.tf_tetra(:,5,ii+5)...
+(ii/domaines.ly)*(-hhy.*domaines.du_tf_tetra(:,5,ii+5)+hhx.*domaines.dv_tf_tetra(:,5,ii+5)));
end;
for ii=-4:4;for jj=-4:4;
a=a+domaines.Bx(ii+5)*domaines.By(jj+5)*hz.* (  ( 4*hhz-2i*pi* ((ii/domaines.lx)*eey-(jj/domaines.ly)*eex )) .*domaines.tf_tetra(:,ii+5,jj+5)...
+hhz.*( (ii/domaines.lx)*domaines.du_tf_tetra(:,ii+5,jj+5)+ (jj/domaines.ly)*domaines.dv_tf_tetra(:,ii+5,jj+5)  )...
-((ii/domaines.lx)*hhx+(jj/domaines.ly)*hhy).*domaines.dw_tf_tetra(:,ii+5,jj+5)) ;
end;end;

end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%                  cylindrique radial                    %%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [a,D]=retcouche_cylindrique_radial(init,w,c);
beta=init{2}(1,:);cao=init{10};c=real(c);
[fmu,fmmu,fep,feep]=retc3_cylindrique_radial(w{1},reshape(w{3},size(w{3},1),size(w{3},3)).',beta);

% retcompare(ab{6}*diag(ab{5}(1:21))*inv(ab{6}),ab1{1}(1:21,1:21)*diag(ab1{5}(1:21))*inv(ab1{1}(1:21,1:21)))
% retcompare(ab{7}*diag(ab{5}(22:41))*inv(ab{7}),ab1{1}(22:41,22:41)*diag(ab1{5}(22:41))*inv(ab1{1}(22:41,22:41)))

if length(cao)>0;
K=apodisecao(cao{1})*diag(beta);% avec apodisation 
else K=diag(beta);end;

ind=rettoeplitz(1:length(fep));
ep=rettoeplitz(fep,ind);Iep=inv(ep);
eep=rettoeplitz(feep,ind);Ieep=inv(eep);
mu=rettoeplitz(fmu,ind);Imu=inv(mu);
mmu=rettoeplitz(fmmu,ind);Immu=inv(mmu);

si=init{8};
if ~isempty(si);si=si(5:8);end;

A_H=calsym((mu-K*Iep*K)*Ieep,si,2,2);[P_H,D_H]=reteig(A_H);
A_E=calsym((ep-K*Imu*K)*Immu,si,0,0);[P_E,D_E]=reteig(A_E);
% pour compatibilité avec le cas classique, on multiplie D_E et D_H par -i Le D du texte theorique est iD du programme  
Z_E=inv(K*Imu*K-ep);Z_H=inv(K*Iep*K-mu);

D_H=retsqrt(-diag(D_H),1);D_E=retsqrt(-diag(D_E),1);D=[D_E;D_H];

Q_EE=.5*calsym(Z_E,si,0,0)*P_E*retdiag(1i*D_E);
Q_EH=.5*calsym(Z_E*K*Imu,si,0,2)*P_H*retdiag(1i*D_H);
Q_HE=.5*calsym(Z_H*K*Iep,si,2,0)*P_E*retdiag(1i*D_E);
Q_HH=.5*calsym(Z_H,si,2,2)*P_H*retdiag(1i*D_H);

if c==1;
N_EE=-.5i*calsym(-Iep+Iep*K*Z_H*K*Iep,si,0,0)*P_E*retdiag(1i*D_E);
N_EH=-.5i*calsym(Iep*K*Z_H,si,0,2)*P_H*retdiag(1i*D_H);
N_HE=-.5i*calsym(Imu*K*Z_E,si,2,0)*P_E*retdiag(1i*D_E);
N_HH=-.5i*calsym(-Imu+Imu*K*Z_E*K*Imu,si,2,2)*P_H*retdiag(1i*D_H);
% O_E=calsym(Ieep,si,0,0)*P_E;% calcul de Dz
% O_H=calsym(Immu,si,2,2)*P_H;% calcul de Bz
O_E=calsym(Immu,si,0,0)*P_E;% calcul de Dz
O_H=calsym(Ieep,si,2,2)*P_H;% calcul de Bz
else;[N_EE,N_EH,N_HE,N_HH,Q_E,Q_H]=deal([]);
end;

parm=struct('dim',init{end}.dim,'genre',2,'type',7,'L',init{end}.L,'nfourier',init{end}.nfourier,'sog',init{end}.sog);
a={Q_HE,Q_HH,Q_EH,Q_EE,D,P_E,P_H,N_EE,N_EH,N_HE,N_HH,w,O_E,O_H,Iep,Imu,Ieep,Immu,K,parm};% w doit etre au numero 12
%    1   2    3    4   5  6   7   8    9   10   11  12 13 14   15  16  17  18  19  end 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [fmu,fmmu,fep,feep]=retc3_cylindrique_radial(x,ep,beta);n=length(beta);alpha=(-n+1:n-1)*(2*pi);
fmu=retf(x,ep(2,:),alpha);fmmu=retf(x,1./ep(1,:),alpha);fep=retf(x,ep(5,:),alpha);feep=retf(x,1./ep(4,:),alpha);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function a=calsym(a,si,n1,n2);
if isempty(si);return;end;
a=si{2+n1}*a*si{1+n2};
