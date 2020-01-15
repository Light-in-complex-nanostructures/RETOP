function [vm,er,ec,ener,enere,enerh,poynting,qvm]=retvm(e,o,xy,wxy,source,k0,champ,parm,fig,texte);
% function [vm,er,ec,ener,enere,enerh,poynting,qvm]=retvm(e,o,xy,wxy,source,k0,champ,parm,fig,texte);
% e,o calcules par retchamp
% xy={x,y}  en 1 D     ou {x,y z} en 2 D
% wxy={wx,wy} en 1 D  ou {wx,wy,wz} en 2 D
% source:position de la 'source' ( max présumé du mode) 
% k0=2*pi/ld dans le cas ou on a un ld complexe
% champ composantes du champ à tracer sur 3 coupes passant par la source
%       parm,fig,texte paramètres de rettchamp pour ce tracé (en 2D création de 3 figures )
% au retour :
% 
% vm(1)=somme(abs(eps*E^2)) /max (abs(eps* abs(E)^2);
% vm(2)=somme(abs( mu*H^2)) /max  (abs(mu* abs(H)^2);
% vm(3)=somme(abs(eps*E^2)+abs(mu*H^2))/max( abs(eps* E^2)+abs(mu* H^2));
% er:  composantes du champ sur la source normalise par le volume reel=somme(abs(eps*E^2)+abs(mu*H^2))
% si il y a  plusieurs points sur la source on prend la moyenne  (valeur principale)
%
% ec:  composantes du champ sur la source normalisé par le volume complexe=somme(eps*E^2+mu*H^2)
% ener: 1/2 * somme(eps*abs(E)^2+mu*abs(H)^2)
% enere: somme(eps*abs(E)^2)
% enerh: somme(mu*abs(H)^2)
%
% poynting flux du vecteur de poynting sortant a travers les faces de la boite
%  poynting(1,1) face x=xmin  poynting(1,2) face x=xmax
%  poynting(2,1) face y=ymin  poynting(2,2) face y=ymax
%  poynting(3,1) face z=zmin  poynting(3,2) face z=zmax
%       ( le maillage doit prendre les valeurs limites (x=retgauss(a,b,n,-m)  etc..)
% qvm:le facteur de qualité déduit
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



if nargin<10;texte='';end;
if nargin<9;fig=0;end;
if nargin<8;parm=0;end;
if nargin<7;champ=[];end;
if nargin<6;k0=1;end;

tchamp=~isempty(champ);

if length(xy)==3;  % 2 D
if nargout>6;poynting=zeros(3,2);calp=1;else;poynting=[];calp=0;end;
x=xy{1};y=xy{2};z=xy{3};wx=wxy{1};wy=wxy{2};wz=wxy{3};
if ~iscell(e);e={e};o={o};z={z};wz={wz};end;
vre=0;vrh=0;vce=0;vch=0;ei=[];ej=[];oi=[];oj=[];zz=[];emax=0;hmax=0;ehmax=0;
% reperage de la source
[prv,isource]=min(abs(x-source(1)));[prv,jsource]=min(abs(y-source(2)));
for ii=1:length(z);zz=[zz,z{ii}];end;[prv,ksource]=min(abs(zz-source(3)));zsource=zz(ksource);
for ii=1:length(e);
oo=retio(o{ii})/k0;ee=retio(e{ii});    
[WX,WZ,WY]=meshgrid(wx,wz{ii},wy);w=WX.*WY.*WZ;clear WX WY WZ;

emax=max(retmax((abs(ee(:,:,:,1)).^2).*abs(oo(:,:,:,4))+...
    (abs(ee(:,:,:,2)).^2).*abs(oo(:,:,:,5))+...
    (abs(ee(:,:,:,3)).^2).*abs(oo(:,:,:,6))),emax);% max de  abs(eps* E^2)

hmax=max(retmax((abs(ee(:,:,:,4)).^2).*abs(oo(:,:,:,1))+...
    (abs(ee(:,:,:,5)).^2).*abs(oo(:,:,:,2))+...
    (abs(ee(:,:,:,6)).^2).*abs(oo(:,:,:,3))),hmax);% max de   abs(mu* H^2)

ehmax=max(retmax((abs(ee(:,:,:,1)).^2).*abs(oo(:,:,:,4))+...
    (abs(ee(:,:,:,2)).^2).*abs(oo(:,:,:,5))+...
    (abs(ee(:,:,:,3)).^2).*abs(oo(:,:,:,6))+...
    (abs(ee(:,:,:,4)).^2).*abs(oo(:,:,:,1))+...
    (abs(ee(:,:,:,5)).^2).*abs(oo(:,:,:,2))+...
    (abs(ee(:,:,:,6)).^2).*abs(oo(:,:,:,3))),ehmax);% max de  abs(eps* E^2)+ abs(mu* H^2)

vre=vre+sum(sum(sum(w.*(abs(ee(:,:,:,1)).^2).*abs(oo(:,:,:,4)))))+...
    sum(sum(sum(w.*(abs(ee(:,:,:,2)).^2).*abs(oo(:,:,:,5)))))+...
    sum(sum(sum(w.*(abs(ee(:,:,:,3)).^2).*abs(oo(:,:,:,6)))));
vrh=vrh+sum(sum(sum(w.*(abs(ee(:,:,:,4)).^2).*abs(oo(:,:,:,1)))))+...
    sum(sum(sum(w.*(abs(ee(:,:,:,5)).^2).*abs(oo(:,:,:,2)))))+...
    sum(sum(sum(w.*(abs(ee(:,:,:,6)).^2).*abs(oo(:,:,:,3)))));
vce=vce+sum(sum(sum(w.*(ee(:,:,:,1).^2).*oo(:,:,:,4))))+...
    sum(sum(sum(w.*(ee(:,:,:,2).^2).*oo(:,:,:,5))))+...
    sum(sum(sum(w.*(ee(:,:,:,3).^2).*oo(:,:,:,6))));
vch=vch+sum(sum(sum(w.*(ee(:,:,:,4).^2).*oo(:,:,:,1))))+...
    sum(sum(sum(w.*(ee(:,:,:,5).^2).*oo(:,:,:,2))))+...
    sum(sum(sum(w.*(ee(:,:,:,6).^2).*oo(:,:,:,3))));
clear w;
if tchamp;
ei=[ei;ee(:,isource,:,:)];ej=[ej;ee(:,:,jsource,:)];oi=[oi;oo(:,isource,:,:)];oj=[oj;oo(:,:,jsource,:)];
end;
k=find(z{ii}==zsource);if ~isempty(k);if tchamp;ek=sum(ee(k,:,:,:),1)/length(k);ok=oo(k,:,:,:);end;es=sum(ee(k,isource(1),jsource(1),:),1)/length(k);end;% source
% si plusieurs points sur la source on prend la moyenne  (valeur principale)


if calp;% flux du vecteur de poynting
poynting(1,1)=poynting(1,1)+retscale(retpoynting(ee(:,1,:,:),[-1,0,0]),(wz{ii}.'*wy));
poynting(1,2)=poynting(1,2)+retscale(retpoynting(ee(:,end,:,:),[1,0,0]),(wz{ii}.'*wy));
poynting(2,1)=poynting(2,1)+retscale(retpoynting(ee(:,:,1,:),[0,-1,0]),(wz{ii}.'*wx));
poynting(2,2)=poynting(2,2)+retscale(retpoynting(ee(:,:,end,:),[0,1,0]),(wz{ii}.'*wx));
if ii==1;poynting(3,1)=retscale(retpoynting(ee(1,:,:,:),[0,0,-1]),(wx.'*wy));end;
if ii==length(e);poynting(3,2)=retscale(retpoynting(ee(end,:,:,:),[0,0,1]),(wx.'*wy));end;
end;

end;  % boucle sur ii

if tchamp;
oi=sqrt(abs(oi(:,:,:,4)));oj=sqrt(abs(oj(:,:,:,4)));ok=sqrt(abs(ok(:,:,:,4)));

end;
if tchamp;% trace du  champ du mode
rettchamp(ei(:,1,:,:)/sqrt(vce-vch),oi,source(1),y,zz,champ,parm,fig,texte);
rettchamp(ej(:,:,1,:)/sqrt(vce-vch),oj,x,source(2),zz,champ,parm,fig+(fig~=0),texte);
rettchamp(ek(1,:,:,:)/sqrt(vce-vch),ok,x,y,source(3),champ,parm,fig+2*(fig~=0),texte);
end;

else; % 1 D
if nargout>6;poynting=zeros(2,2);calp=1;else;poynting=[];calp=0;end;
x=xy{1};y=xy{2};wx=wxy{1};wy=wxy{2};
vre=0;vrh=0;vce=0;vch=0;emax=0;hmax=0;
o=o/k0;
% reperage de la source
[prv,isource]=min(abs(x-source(1)));[prv,jsource]=min(abs(y-source(2)));ysource=y(jsource);
es=e(jsource,isource,:);
emax=retmax((abs(e(:,:,1)).^2).*abs(o(:,:,1))); % max abs(eps* E^2)
hmax=retmax((abs(e(:,:,2)).^2).*abs(o(:,:,2))+(abs(e(:,:,3)).^2).*abs(o(:,:,3))); % max abs(mu* H^2)
ehmax=retmax((abs(e(:,:,1)).^2).*abs(o(:,:,1))+(abs(e(:,:,2)).^2).*abs(o(:,:,2))+(abs(e(:,:,3)).^2).*abs(o(:,:,3)));% max (abs(eps* E^2) + abs(mu* H^2))
w=wy.'*wx;
vre=sum(sum(w.*(abs(e(:,:,1)).^2).*abs(o(:,:,1))));
vrh=sum(sum(w.*(abs(e(:,:,2)).^2).*abs(o(:,:,2))))+sum(sum(w.*(abs(e(:,:,3)).^2).*abs(o(:,:,3))));
vce=sum(sum(w.*(e(:,:,1)).^2.*o(:,:,1)));
vch=sum(sum(w.*(e(:,:,2)).^2.*o(:,:,2)))+sum(sum(w.*(e(:,:,3)).^2.*o(:,:,3)));
clear w;
o=sqrt(abs(o(:,:,1).*o(:,:,2)));

if calp;% flux du vecteur de poynting
poynting(1,1)=sum(retpoynting(e(:,1,:),[-1,0]).*wy.');
poynting(1,2)=sum(retpoynting(e(:,end,:),[1,0]).*wy.');
poynting(2,1)=sum(retpoynting(e(1,:,:),[0,-1]).*wx);
poynting(2,2)=sum(retpoynting(e(end,:,:),[0,1]).*wx);
end;




if tchamp;rettchamp(e/sqrt(vce-vch),o,x,y,0,champ,0,fig,texte);end;% trace du champ du mode
end;  % fin 1 D

ener=(vre+vrh)/2;enere=vre;enerh=vrh;
er=i*squeeze(es/sqrt(vre+vrh));% champ normalise par volume modal reel
ec=squeeze(es/sqrt(vce-vch));% champ normalise par volume modal complexe
vm=zeros(1,3);
vm(3)=(vre+vrh)/ehmax;vm(1)=vre/emax;vm(2)=vrh/hmax;



if calp;ldz=2*pi/(real(k0)-i*sum(sum(poynting))/enerh);qvm=real(ldz)/(2*imag(ldz));end; % facteur de qualite

