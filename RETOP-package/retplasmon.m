function varargout=retplasmon(varargin);

%
%
%          %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%          %  PARAMETRES D'UN PLASMON  %
%          %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% function  plasmon=retplasmon(indice_metal,indice_dielectrique,lambda,betay);
% retourne une structure contenant;
%      constante_de_propagation, khi_metal, khi_dielectrique
%      test_d_existence = 1 si le plasmon existe ,0 sinon
% si lambda est precise on a aussi les distances d'attenuation en intensite d'un facteur e
%      distance_de_propagation,penetration_metal,penetration_dielectrique
%      attenuation en db/mm si lambda est en microns
%      poynting=flux du vecteur de Poynting d'un plasmon dont le vecteur H vaut 1 sur le dioptre
%         (poynting = poynting_dielectrique + poynting_air )
%      int_lorentz: integrale de Lorentz de ce meme plasmon ( poynting voisin de int_lorentz/4 )
%      res: residu( (nm*nd)^3/(nd^4-nm^4) (res=i/int_lorentz)
%
% par defaut   indice_dielectrique=1
%
%% EXEMPLE  distance de propagation du plasmon de l'or fonction de lambda
%  figure;ld=linspace(.4,2.5,1000);
%  pl=retplasmon(retindice(ld,2),1,ld);plot(ld,pl.distance_de_propagation);
%  grid;xlabel('lambda microns ');title('distance de propagation du plasmon de l''or (microns)')
%
%
%
%          %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%          %  PARAMETRES D'UNE PARTICULE  %
%          %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% particule=retplasmon(L,H,k0,n_dielectrique,n_metal,n_particule,teta,m)
%
%                si H fini   |                                |     si H =inf 
%    si H>0                  |          si H<0                |
%                            |             L                  |
%            n_dielectrique  |            XXXX n_particule    |
%                            |            XXXX -H             |
%          L                 |            XXXX                |            L
% *******    *********       |    *********************       |    *******    *********
% *******  H *********       |    *********************       |    *******    *********
% ********************       |    ***   n_metal  *****        |    *******    *********
% ********************       |    *********************       |    *******    *********
% Entrees: k0=2*pi/ld  teta en radians  m nombre de termes de fourier
%   si m a une partie complexe le nombre de termes de Fourier est real(m),
%   et on ne calcule pas la reponse de la particule au plasmon (le couplage plasmon --> mode est donné par la reciprocité)
%
% Sortie:  particule est une structure à champs
% 
% particule = 
%              teta_degres: 0
%                       LF: 0.3200
%                       HF: Inf
%                       LD: 0.8000
%                  n_metal: 0.1803 + 5.1287i
%           n_dielectrique: 1
%              n_particule: 1
%             poynting_inc: 0.1315
%            poynting_mode: 0.1892
%         poynting_plasmon: 0.1632
%            Mode_incident: [1x1 struct] (seulement si H=inf)
%     Onde_plane_incidente: [1x1 struct]
% 
% ans = 
% 						amplitude_mode_incident: 1
% 							  amplitude_plg: 0.3581 - 0.0131i
% 							  amplitude_pld: 0.3581 - 0.0131i
% 							 amplitude_mode_reflechi: 0.0255 - 0.2329i
%           Mode_incident=  Energie_Mode_Incident: 0.1892
% 								Energie_Plg: 0.0210
% 								Energie_Pld: 0.0210
% 							   Energie_Mode_Reflechi: 0.0104
% 	
%                                                  amplitude_incidente: 0.5034 - 0.0974i
%                                                        amplitude_plg: -0.1713 - 0.0305i
%                                                        amplitude_pld: -0.1713 - 0.0305i
% Onde_plane_incidente=                                 amplitude_mode: 0.4446 + 0.0020i
%                         Energie_Incidente_sur_la_largeur_de_la_particule: 0.0421
%                                                          Energie_Plg: 0.0049
%                                                          Energie_Pld: 0.0049
%                                                         Energie_Mode: 0.0374
%
%% EXEMPLE
%% m=50;teta=0;LD=.6328;k0=2*pi/LD;L=.4*LD;n_metal=retindice(LD,2);n_dielectrique=1;n_particule=1;H=inf;
% m=150;teta=0;LD=.6328;k0=2*pi/LD;L=.25*LD;n_metal=retindice(LD,2);n_dielectrique=1;n_particule=1;H=inf;
% particule=retplasmon(L,H,k0,n_dielectrique,n_metal,n_particule,teta,m)
% disp('particule.Onde_plane_incidente='); disp(particule.Onde_plane_incidente)
% if ~isfinite(H);disp('particule.Mode_incident='); disp(particule.Mode_incident),end;
% disp('particule.Plasmon_incident='); disp(particule.Plasmon_incident)
% Beta_2=particule.Onde_plane_incidente.Energie_Plg/particule.Onde_plane_incidente.Energie_Incidente_sur_la_largeur_de_la_particule
% if ~isfinite(H);
% Alpha_mode_plasmon_2=particule.Mode_incident.Energie_Plg/particule.Mode_incident.Energie_Mode_Incident
% Alpha_plasmon_mode_2=particule.Plasmon_incident.Energie_Mode/particule.Plasmon_incident.Energie_Plasmon_Incident
% T_2=particule.Onde_plane_incidente.Energie_Mode/particule.Onde_plane_incidente.Energie_Incidente_sur_la_largeur_de_la_particule
% end
%
% See also:RETPL,RETDB,RETPARTICULE


if nargin<=4;[varargout{1:nargout}]=plasmon(varargin{:});
else;[varargout{1:nargout}]=cal_particule(varargin{:});
end;                                                    

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function  pl=plasmon(nm,nd,ld,betay);
st_warning=warning;warning off;
if nargin<4;betay=0;end;
if nargin<2;nd=1;end; % air par defaut
beta=sqrt(1./(nd.^-2+nm.^-2)-betay^2);khim=retsqrt(nm.^2-(beta.^2+betay^2),-1);khid=retsqrt(nd.^2-(beta.^2+betay^2),-1);
test=abs((khim+(nm./nd).^2.*khid)./(2*khim))<.1;
%e=[beta,betay,-i*(beta^2+betay^2)/khim];
pl=struct('constante_de_propagation',beta,'khi_metal',khim,'khi_dielectrique',khid,'test_d_existence',test);
if nargin>2;a=ld/(4*pi);
pl.distance_de_propagation=a./abs(imag(beta));
pl.penetration_metal=a./abs(imag(khim));
pl.penetration_dielectrique=a./abs(imag(khid));
pl.db_par_mm=retdb(beta,ld);
k0=2*pi./ld;
pl.poynting=real((beta.^2+betay^2)./(4*k0.*nm.^2.*beta))./imag(khim)+real((beta.^2+betay^2)./(4*k0.*nd.^2.*beta))./imag(khid);
pl.poynting_metal=real((beta.^2+betay^2)./(4*k0.*nm.^2.*beta))./imag(khim);
pl.poynting_dielectrique=real((beta.^2+betay^2)./(4*k0.*nd.^2.*beta))./imag(khid);
pl.int_lorentz=i*((beta.^2+betay^2)./beta).*(1./(k0.*nm.^2.*khim)+1./(k0.*nd.^2.*khid));
pl.res=(nm.*nd).^3./(nd.^4-nm.^4);
end;
warning(st_warning);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function particule=cal_particule(L,H,k0,n_dielectrique,n_metal,n_particule,teta,m);
X=0;
LD=2*pi/k0;
beta=k0*n_dielectrique*sin(teta);
plasmon=retplasmon(n_metal,n_dielectrique,LD);
r=retabeles(2,[n_dielectrique,n_metal],0,beta/k0);
mode=retpmode(2,n_metal,n_dielectrique,L,k0);

[S,Vp,poynting_inc,poynting_mode]=retparticule_bis(X,L,H,k0,n_dielectrique,n_metal,n_particule,beta,real(m),1,inf);
Dif=reteval(S);
[plg_mode,pld_mode,Plg_mode,Pld_mode]=retparticule_bis(n_dielectrique,n_metal,X,Dif(:,1),[-inf,inf],k0);  % amplitude sur les plasmons mode incident dans la fente
[plg_op,pld_op,Plg_op,Pld_op]=retparticule_bis(n_dielectrique,n_metal,X,Dif(:,2),[-inf,inf],k0);% amplitude sur les plasmons onde plane incidente

if imag(m)==0;
beta_pl=k0*plasmon.constante_de_propagation;
S_pl=retparticule_bis(X,L,H,k0,n_dielectrique,n_metal,n_particule,beta_pl,real(m),1,inf);
Dif_pl=reteval(S_pl);
[plg_pl,pld_pl,Plg_pl,Pld_pl]=retparticule_bis(n_dielectrique,n_metal,X,Dif_pl(:,2),[-inf,inf],k0);  % amplitude sur les plasmons plasmon incident
% [Ep,Em,Angles]=retparticule_bis(n_dielectrique,n_metal,Dif_pl,u,X,k0,1);
else;
Dif_pl=nan(5,2);[plg_pl,pld_pl,Plg_pl,Pld_pl]=deal(nan(1,2));
end;



particule=struct('teta_degres',teta*180/pi,'L',L,'H',H,'LD',LD,'n_metal',n_metal,'n_dielectrique',n_dielectrique,'n_particule',n_particule,'poynting_inc',poynting_inc,'poynting_plasmon',plasmon.poynting,'int_lorentz_plasmon',plasmon.int_lorentz);
if isfinite(H);
particule.Onde_plane_incidente=struct('amplitude_incidente',1/(1+r),'amplitude_plg',plg_op(1),'amplitude_pld',pld_op(2),'Energie_Incidente_sur_la_largeur_de_la_particule',poynting_inc*L,'Energie_Plg',Plg_op(1),'Energie_Pld',Pld_op(2));
if imag(m)==0;particule.Plasmon_incident=struct('amplitude_pl_incident',1,'amplitude_plasmon_reflechi',plg_pl(1),'amplitude_plasmon_transmis',1+pld_pl(2),'Energie_Plasmon_Incident',plasmon.poynting,'Energie_Plasmon_Reflechi',plasmon.poynting*Plg_pl(1),'Energie_Plasmon_Transmis',plasmon.poynting*abs(1+pld_pl(2))^2);end;
else;
mode=retpmode(2,n_metal,n_dielectrique,L,k0);particule.int_lorentz_mode=mode.int_lorentz;particule.test_reciprocite=nan;
particule.Mode_incident=struct('amplitude_mode_incident',1,'amplitude_plg',plg_mode(1),'amplitude_pld',pld_mode(2),'amplitude_mode_reflechi',Dif(5,1),'Energie_Mode_Incident',poynting_mode,'Energie_Plg',Plg_mode(1),'Energie_Pld',Pld_mode(2),'Energie_Mode_Reflechi',poynting_mode*abs(Dif(5,1))^2);
particule.Onde_plane_incidente=struct('amplitude_incidente',1/(1+r),'amplitude_plg',plg_op(1),'amplitude_pld',pld_op(2),'amplitude_mode',Dif(5,2),'Energie_Incidente_sur_la_largeur_de_la_particule',poynting_inc*L,'Energie_Plg',Plg_op(1),'Energie_Pld',Pld_op(2),'Energie_Mode',poynting_mode*abs(Dif(5,2))^2);
if imag(m)==0;
particule.Plasmon_incident=struct('amplitude_pl_incident',1,'amplitude_plasmon_reflechi',plg_pl(1),'amplitude_plasmon_transmis',1+pld_pl(2),'amplitude_mode',Dif_pl(5,2),'Energie_Plasmon_Incident',plasmon.poynting,'Energie_Plasmon_Reflechi',plasmon.poynting*Plg_pl(1),'Energie_Plasmon_Transmis',plasmon.poynting*abs(1+pld_pl(2))^2,'Energie_Mode',poynting_mode*abs(Dif_pl(5,2))^2);
particule.test_reciprocite=abs((particule.Mode_incident.amplitude_plg*particule.int_lorentz_plasmon)/(particule.Plasmon_incident.amplitude_mode*particule.int_lorentz_mode)-1);
else;
prv=particule.Mode_incident.amplitude_plg*(particule.int_lorentz_plasmon/particule.int_lorentz_mode);	
particule.Plasmon_incident=struct('amplitude_pl_incident',1,'amplitude_mode',prv,'Energie_Plasmon_Incident',plasmon.poynting,'Energie_Mode',poynting_mode*abs(prv));
end;
end;
