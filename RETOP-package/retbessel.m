function varargout=retbessel(type,nu,varargin);
% function B=retbessel(type,mêmes variables que dans Matlab);
% Accélération du calcul des fonctions de bessel d'indice entier<0
%
%  Curieusement le calcul des fonctions de bessel d'indice entier <0 avec Matlab ( 7.3 )
%  est beaucoup plus long (facteur 9) que celui des fonctions d'indice entier >0
% En utilisant la relation: (sauf pour les fonctions de Kelvin) bessel(-nu,x)=(-1)^nu   bessel(nu,x)  
% on obtient un gain de temps pour les nu<0
% Remarque: bien que la propriete ci dessus ne soit vraie que pour les indices entiers,
% le programme est prevu pour fonctionner pour les indices non entiers 
% type= 'j','y','k','h'
%
% ATTENTION: nu et la variable Z doivent etre respectivement une ligne et une colonne
% Le résultat est de dimension [length(Z),length(nu)] 
% 
%% EXEMPLE
% L=-20:20;x=randn(5000,1)+i*randn(5000,1);[xx,LL]=ndgrid(x,L);
% tic;H=besselh(LL,1,xx,1);cpu_bessel=toc
% tic;HH=retbessel('h',L,1,x,1);cpu_retbessel=toc, retcompare(H,HH)
%
% See also: BESSELJ BESSELY BESSELH BESSELK

% if all(type=='h');nx=2;else;nx=1;end;
% sz_x=size(varargin{nx});num_x=prod(sz_x);
% sz_nu=size(nu);num_nu=prod(sz_nu);
% if num_x==1 & num_nu>1;varargin{nx}=varargin{nx}*ones(sz_nu);end
% if num_x>1 & num_nu==1;nu=nu*ones(sz_x);end
% if (sz_x(1)==1 & sz_x(2)>1) & (sz_nu(2)==1 & sz_nu(1)>1) | (sz_x(2)==1 & sz_x(1)>1) & (sz_nu(1)==1 & sz_nu(2)>1) ;[varargin{nx},nu]=ndgrid(varargin{nx},nu);end


g=find(mod(nu,1)==0);

ff=g(nu(g)<0);
% if isempty(ff);mu=nu;
% else;
fff=ff(mod(-nu(ff),2)==1);
mu=nu;mu(ff)=abs(mu(ff));
[mu,k,kk]=retelimine(mu,-1);% tri des valeurs à calculer
%end;

if retversion<8;
%if vv<=6;
switch(type);
case 'j';[varargout{1:nargout}]=besselj(mu,varargin{:});	
case 'y';[varargout{1:nargout}]=bessely(mu,varargin{:});	
case 'k';[varargout{1:nargout}]=besselk(mu,varargin{:});	
case 'h';[varargout{1:nargout}]=besselh(mu,varargin{:});	
end
%if isempty(ff);return;end;

    else;% bugg 7 2013 Anthony
    if nargout>1;varargout{2}=0;end    
    %if type=='h';
    if length(mu)>1;
    if type=='h';num=2;else;num=1;end;
    if length(varargin{num})~=1;[varargin{num},mu]=ndgrid(varargin{num},mu);end;
    end;
    switch(type);
    case 'j';varargout{1}=besselj(mu,varargin{:});	
    case 'y';varargout{1}=bessely(mu,varargin{:});	
    case 'k';varargout{1}=besselk(mu,varargin{:});	
    case 'h';varargout{1}=besselh(mu,varargin{:});	
    end
    end;
if ~isempty(varargout{1});
varargout{1}=varargout{1}(:,kk);
if type~='k';varargout{1}(:,fff)=-varargout{1}(:,fff);end;
end;









