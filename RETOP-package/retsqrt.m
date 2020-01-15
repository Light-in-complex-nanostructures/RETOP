function [y,z]=retsqrt(x,k,tet,ww);
%  function y=retsqrt(x,k,tet,w);
%  determination de sqrt(x)
%
%  k=0 ou abs: DETERMINATION DES ONDES PLANES   PRIORITEE A L'ATTENUATION   (coupure de Petit) 
%     imag(y)>0 et si imag(y)==0,real(y)>=0  
%     exp(i y *z) onde bornee si z>0 ou  dirigee vers les z>0
%
%  k=-1 DETERMINATION DES ONDES PLANES   exp( iy *z ) est une onde dirigee vers les  z>0
%  tet :angle de la coupure (rad) pi/2 par defaut (coupure de Maystre )
%
%  k=1 DETERMINATION DES MODES   exp( y *z ) est une onde dirigee vers les  z<0
%    ( si yy=iy      exp(i yy *z) est une onde dirigee vers les  z>0  )
%  tet :angle de la coupure (rad) pi/2 par defaut (coupure de Maystre )
%
% Majoration de Wojciech pour les cas trop amplifiés
% si imag(y)<-w (pour k=0 ou -1) ou  real(y)>w (pour k=1), on change le signe de y
% par defaut w=inf
%
%  pour changer tet  retsqrt(0,0,tet)  tet reste alors  changé pour la suite
%  pour changer tet et w  retsqrt(0,0,tet,w)  tet et restent alors  changé pour la suite (tant que l'on ne fait pas 'ret' ou 'clear retsqrt')
%   la borne w doit etre adaptee aux unités: c'est une borne sur les valeurs propres donc k0 neff
%
%  [tet,w]=retsqrt retourne la valeur actuelle de tet et w
%
%% Exemples et test
% [X,Y]=ndgrid(linspace(-2,2,51),linspace(-2,2,53));Z=X+i*Y;figure;
% subplot(2,2,1);plot(retsqrt(Z),'.k');title('retsqrt(Z)');axis equal
% subplot(2,2,2);plot(retsqrt(Z,-1),'.k');title('retsqrt(Z,-1)');axis equal
% subplot(2,2,3);plot(retsqrt(Z,1),'.k');title('retsqrt(Z,1)');axis equal
% subplot(2,2,4);retsqrt(0,0,pi/2,.5);plot(retsqrt(Z,1),'.k');title('retsqrt(0,0,pi/2,.5);retsqrt(Z,1)');axis equal;retsqrt(0,0,pi/2,inf);
%
% See also: RETTNEFF


persistent teta e1 e2 w;if isempty(teta);teta=pi/2;e1=-i;e2=(1-i)/rettsqrt(2);w=inf;end;
if nargin>2;% modification des parametres par defaut
if nargin<4;ww=inf;end;
teta=tet;w=ww;
if teta~=pi/2;e1=exp(-i*teta);e2=exp(-.5i*teta);else;e1=-i;e2=(1-i)/rettsqrt(2);end;
end;
if nargin<1;y=teta;z=w;return;end;
if nargin<2;k=0;end;
switch(k);
case 0;   % coupure de Petit
y=i*conj(rettsqrt(-conj(x)));
if isfinite(w);f=find(imag(y)<-w);y(f)=-y(f);end;% Majoration de Wojciech
case 1;  % coupure a teta determination des MODES
if teta==pi;y=-i*rettsqrt(-x);else;y=conj(rettsqrt(e1*conj(x)))*e2;end;
if isfinite(w);f=find(real(y)<-w);y(f)=-y(f);end;% Majoration de Wojciech
case -1;  % coupure a teta determination des ONDES PLANES
if teta==pi;y=rettsqrt(x);else;y=i*conj(rettsqrt(-e1*conj(x)))*e2;end;
if isfinite(w);f=find(imag(y)<-w);y(f)=-y(f);end;% Majoration de Wojciech
end;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function x=rettsqrt(x);% c'est la racine mathematique
if isempty(x);return;end;
x=sqrt(x);f=find((real(x)<0)|(real(x)==0&imag(x)<0));x(f)=-x(f);