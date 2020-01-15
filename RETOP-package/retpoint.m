function varargout=retpoint(varargin);

% e=retpoint(ep,mu,s,centre,x,y,z,clone);
% champ cree par une source ponctuelle s au point: centre ,dans un milieu isotrope ep mu
%   en 2 D          e(z,x,y,1:6)  e=retpoint(ep,mu,s,centre,x,y,z,clone);      <-----------length(centre)=3
%   en 1 D conique  e(z,x,y,1:6)    e=retpoint(ep,mu,gama,s,centre,x,y,z,clone); <-----------length(centre)=3,length(S)=6
%   en 1 D          e(y,x,1:3)    e=retpoint(ep,mu,s,centre,x,y,pol,clone);    <-----------length(centre)=2
%   en 0 D          e(x,1:2)      e=retpoint(ep,mu,s,centre,x,pol,clone);      <-----------length(centre)=1
%
% quand pol=2 ep et mu sont ÈchangÈs dans le programme, il faut mettre en entree ep et mu physiques: ep=k0*n^2 ,mu=k0 
% la definition de e est conforme a celle de retchamp
%   dans le cas de polarisation E// (pol=0)  ez=e(:,:,1),hx=e(:,:,2),hy=e(:,:,3) 
%   dans le cas de polarisation H// (pol=2)  hz=e(:,:,1),ex=-e(:,:,2),ey=-e(:,:,3) 
% en 2D: ex=e(:,:,:,1),ey=e(:,:,:,2),ez=e(:,:,:,3),hx=e(:,:,:,4),hy=e(:,:,:,5),hz=e(:,:,:,6)
%
%   clone: 1 clonage 0 pas clonage(0 par defaut)
%   si clone a une partie imaginaire x y z ont mÍme dimension et ne seront pas 'meshÈs'
%            (alors [Z,X,Y]=ndgrid(z,x,y) en 2D [Y,X]=ndgrid(y,x) en 1D )
%
% Sources pseudo pÈriodiques: utilise une tabulation du champ calculÈ par FFT sur une hauteur zmax 
% e=retpoint(ep,mu,s,centre,x,y,z,  d,beta0,zmax,clone); d periode, beta0 incidence
%    e est le champ total (on ajoute la source centrale comme dans retbahar)
% (si on veut le champ moins la source, il faut ajouter le parametre parm=struct('point',0)) ) 
%
% Il est possible de traiter simultanement plusieurs sources
% size(S)=[3,nb_Sources];en 1D   size(S)=[6,nb_Sources];en 1D
% alors e a une dimension de plus (nb_Sources)
%
% % %  Exemple 3D (rÈseau 2D)
% ld=.8;  k0=2*pi/ld; n=1.5;S=[1,0,0,0,0,0];  centre=[0,0,0]  ;
% x=linspace(-2*ld,2*ld,101);y=linspace(-2*ld,2*ld,105);  z=ld/4; 
% e=retpoint(k0*n^2,k0,S, centre,x,y,z);rettchamp(e,[],x,y,z,[1,14]);
%
% 
% %  Exemple 2D (rÈseau 1D)
% pol=2;ld=.8;  k0=2*pi/ld; n=1.5; S=[0,1,0];  centre=[0,0]  ;
% x=linspace(-2*ld,2*ld,101);y=linspace(-2*ld,2*ld,105);  
% e=retpoint(k0*n^2,k0,S, centre,x,y,pol);rettchamp(e,[],x,y,pol,[1,14]);%
% 
% %  Exemple 1D (rÈseau 0D)
% pol=2;ld=.8;  k0=2*pi/ld; n=1.5; S=[0,1];  centre=0  ;
% x=linspace(-2*ld,2*ld,101) ;
% e=retpoint(k0*n^2,k0,S, centre,x,pol);figure;plot(x,real(e(:,1)));
%
% % 
% %%   Exemple: energie emise dans le bulk:
%   S3=randn(6,1)+i*randn(6,1); S2=randn(3,1)+i*randn(3,1);S3_cl=S3.*[1i;1i;1i;1;1;1];
%   S1=randn(2,1)+i*randn(2,1);
%   k0=2*pi/1.5;n=3.5;gama=.7*k0*n;
% %% S2 et S3 sont considerees comme non clonÈes, S3_cl est clonÈ
% %% pour 1D 2D et 3D il n'intervient que les valeurs absolues, pour le conique il y a un terme de couplage E H
%  %% En clonant
%   Ener_3D   =.5*imag(retpoint(k0*n^2,k0,S3,[0,0,0],0,0,0,1+i)*conj(S3));
%   Ener_E_2D =.5*imag(retpoint(k0*n^2,k0,S2,[0,0],0,0,0,1+i)*conj(S2));
%   Ener_H_2D =.5*imag(retpoint(k0*n^2,k0,S2,[0,0],0,0,2,1+i)*conj(S2));
%   Ener_E_1D =.5*imag(retpoint(k0*n^2,k0,S1,0,0,0,1+i)*conj(S1));
%   Ener_H_1D =.5*imag(retpoint(k0*n^2,k0,S1,0,0,2,1+i)*conj(S1));
%   Ener_conique=.5*imag(retpoint(k0*n^2,k0,gama,S3_cl,[0,0,0],0,0,0,1+i)*conj(S3_cl));
%  
%   disp(rettexte(Ener_3D,Ener_E_2D,Ener_H_2D,Ener_E_1D,Ener_H_1D,Ener_conique));
%  %% En ne clonant pas
%   Ener_3D   =.5*real(retpoint(k0*n^2,k0,S3,[0,0,0],0,0,0,i)*conj(S3.*[-1;-1;-1;1;1;1]));
%   Ener_E_2D =.5*real(retpoint(k0*n^2,k0,S2,[0,0],0,0,0,i)*conj(S2.*[-1;1;1]));
%   Ener_H_2D =.5*real(retpoint(k0*n^2,k0,S2,[0,0],0,0,2,i)*conj(S2.*[-1;1;1]));
%   Ener_E_1D =.5*real(retpoint(k0*n^2,k0,S1,0,0,0,i)*conj(S1.*[-1;1]));
%   Ener_H_1D =.5*real(retpoint(k0*n^2,k0,S1,0,0,2,i)*conj(S1.*[-1;1]));
%   Ener_conique=.5*real(retpoint(k0*n^2,k0,gama,S3,[0,0,0],0,0,0,i)*conj(S3.*[-1;-1;-1;1;1;1]));
%  
%   disp(rettexte(Ener_3D,Ener_E_2D,Ener_H_2D,Ener_E_1D,Ener_H_1D,Ener_conique));
%  %% Formules thÈoriques
%   Ener_3D   = k0^2*n/(12*pi)*sum(abs(S3(1:3)).^2)+k0^2*n^3/(12*pi)*sum(abs(S3(4:6)).^2);
%   Ener_E_2D =(k0/8)*abs(S2(1)).^2+(k0/16)*n^2*sum(abs(S2(2:3)).^2);
%   Ener_H_2D =(k0/8)*n^2*abs(S2(1)).^2+(k0/16)*sum(abs(S2(2:3)).^2);
%   Ener_E_1D =1/(4*n)*abs(S1(1))^2+n/4*abs(S1(2))^2;
%   Ener_H_1D =n/4*abs(S1(1))^2+1/(4*n)*abs(S1(2))^2;
%   Ener_conique=k0*((n^2-(gama/k0)^2)/8)*(abs(S3(3))^2/n^2+abs(S3(6))^2)+k0*((n^2+(gama/k0)^2)/16)*(abs(S3(1)^2)/n^2+abs(S3(2)^2)/n^2+abs(S3(4)^2)+abs(S3(5)^2))-(gama/4)*real(conj(S3(1))*S3(5)-conj(S3(2))*S3(4));
%   Ener_conique_cl=k0*((n^2-(gama/k0)^2)/8)*(abs(S3_cl(3))^2/n^2+abs(S3_cl(6))^2)+k0*((n^2+(gama/k0)^2)/16)*(abs(S3_cl(1)^2)/n^2+abs(S3_cl(2)^2)/n^2+abs(S3_cl(4)^2)+abs(S3_cl(5)^2))+(gama/4)*imag(conj(S3_cl(1))*S3_cl(5)-conj(S3_cl(2))*S3_cl(4));
%   disp(rettexte(Ener_3D,Ener_E_2D,Ener_H_2D,Ener_E_1D,Ener_H_1D,Ener_conique,Ener_conique_cl));
%  %
%%% See also:RETCHAMP RETSOMMERFELD RETS RETBAHAR

% retpoint(n^2,1,S,k0*centre,k0*x,k0*y,k0*z)=retpoint(k0*n^2,k0,S/k0^2,centre,x,y,z); % 2D
% retpoint(n^2,1,S,k0*centre,k0*x,k0*y,pol)=retpoint(k0*n^2,k0,S/k0,centre,x,y,pol);  % 1D


%%%%%%%%%%%%%%%%%%%%
%   AIGUILLAGE     %
%%%%%%%%%%%%%%%%%%%%


if length(varargin{4})>=6;
%if nargin<10;[varargout{1:nargout}]=retpoint_conique(varargin{:});else;[varargout{1:nargout}]=retpoint_conique(varargin{:});end;	%??
[varargout{1:nargout}]=retpoint_conique(varargin{:});	
else

if nargin>8 ;
if length(varargin{4})==3;[varargout{1:nargout}]=retpoint_periodique_2D(varargin{:});else;[varargout{1:nargout}]=retpoint_periodique_1D(varargin{:});end;
else;
[varargout{1:nargout}]=retpoint_012D(varargin{:});
end;
end;

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function e=retpoint_conique1(ep,mu,gama,s,xs,x,y,z,clone);
% khi2=ep*mu-gama^2;khi=retsqrt(khi2,-1);
% if nargin<9;clone=0;end;xmesh=imag(clone)~=0;clone=real(clone);
% if clone==0;s(1:3)=1i*s(1:3);end;
% if xmesh;xx=x;yy=y;zz=z;else;[zz,xx,yy]=ndgrid(z,x,y);end;
% xx=xx(:)-xs(1);yy=yy(:)-xs(2);zz=zz(:)-xs(3);
% r=sqrt(xx.^2+yy.^2);
% 
% e=zeros(numel(r),6);
% [h0,h1sr,h2sr2]=calh_conique(khi*r);
% e(:,3)=(-.25i*khi2)*( (-s(3)/ep)*h0+(xx*s(5)-yy*s(4)+(1i*gama/ep)*(xx*s(1)+yy*s(2))).*h1sr);
% e(:,6)=(-.25i*khi2)*( (-s(6)/mu)*h0+(xx*s(2)-yy*s(1)+(1i*gama/mu)*(xx*s(4)+yy*s(5))).*h1sr);
% e(:,1)=(-.25i*mu)*( (yy*((khi2/mu)*s(6))-s(1)+(1i*gama/mu)*s(5)).*h1sr -khi2*(xx.*yy*s(2)-yy.^2*s(1)+(1i*gama/mu)*(xx.*yy*s(4)+yy.^2*s(5))  ).*h2sr2)...
% 	+(gama/4)* ( (xx*((khi2/ep)*s(3))+s(5)+(1i*gama/ep)*s(1)).*h1sr -khi2*(xx.^2*s(5)-xx.*yy*s(4)+(1i*gama/ep)*(xx.^2*s(1)+xx.*yy*s(2))  ).*h2sr2);
% e(:,2)=(gama/4)* ( (yy*((khi2/ep)*s(3))-s(4)+(1i*gama/ep)*s(2)).*h1sr -khi2*(xx.*yy*s(5)-yy.^2*s(4)+(1i*gama/ep)*(xx.*yy*s(1)+yy.^2*s(2))  ).*h2sr2)...
%     +(.25i*mu)*( (xx*((khi2/mu)*s(6))+s(2)+(1i*gama/mu)*s(4)).*h1sr -khi2*(xx.^2.*s(2)-xx.*yy*s(1)+(1i*gama/mu)*(xx.^2*s(4)+xx.*yy*s(5))  ).*h2sr2);
% 
% 
% e(:,4)=(-.25i*ep)*( (yy*((khi2/ep)*s(3))-s(4)+(1i*gama/ep)*s(2)).*h1sr -khi2*(   xx.*yy*s(5)-yy.^2*s(4)+(1i*gama/ep)*(xx.*yy*s(1)+yy.^2*s(2))  ).*h2sr2)...
% 	+(gama/4)* ( (xx*((khi2/mu)*s(6))+s(2)+(1i*gama/mu)*s(4)).*h1sr -khi2*(xx.^2*s(2)-xx.*yy*s(1)+(1i*gama/mu)*(xx.^2*s(4)+xx.*yy*s(5))  ).*h2sr2);
% e(:,5)=(gama/4)* ( (yy*((khi2/mu)*s(6))-s(1)+(1i*gama/mu)*s(5)).*h1sr -khi2*(xx.*yy*s(2)-yy.^2*s(1)+(1i*gama/mu)*(xx.*yy*s(4)+yy.^2*s(5))  ).*h2sr2)...
%     +(.25i*ep)*( (xx*((khi2/ep)*s(3))+s(5)+(1i*gama/ep)*s(1)).*h1sr -khi2*(xx.^2.*s(5)-xx.*yy*s(4)+(1i*gama/ep)*(xx.^2*s(1)+xx.*yy*s(2))  ).*h2sr2);
% prv=exp(1i*gama*zz);for ii=1:6;e(:,ii)=e(:,ii).*prv;end;
% if khi2<0;e(r==0,:)=0;end;% energie sur la source
% if clone==0;e(:,4:6)=-1i*e(:,4:6);end;
% if ~xmesh;e=reshape(e,length(z),length(x),length(y),6);end;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function e=retpoint_conique(ep,mu,gama,s,xs,x,y,z,clone);
% Il est possible de traiter simultanÈment plusieurs sources
% size(S)=[3,nb_Sources];en 1D   size(S)=[6,nb_Sources];en 1D
% alors e a une dimension de plus (nb_Sources)
if size(s,1)==1 | size(s,2)==1 ;s=s(:);end; 

khi2=ep*mu-gama^2;khi=retsqrt(khi2,-1);
if nargin<9;clone=0;end;xmesh=imag(clone)~=0;clone=real(clone);
if clone==0;s(1:3,:)=1i*s(1:3,:);end;
if xmesh;xx=x;yy=y;zz=z;else;[zz,xx,yy]=ndgrid(z,x,y);end;
xx=xx(:)-xs(1);yy=yy(:)-xs(2);zz=zz(:)-xs(3);
r=sqrt(xx.^2+yy.^2);

e=zeros(numel(r),size(s,2),6);
[h0,h1sr,h2sr2]=calh_conique(khi*r);
prv=exp(1i*gama*zz);h0=prv.*h0;h1sr=prv.*h1sr;h2sr2=prv.*h2sr2;
xxh1sr=xx.*h1sr;yyh1sr=yy.*h1sr;
xx2h2sr2=xx.^2.*h2sr2;yy2h2sr2=yy.^2.*h2sr2;xxyyh2sr2=xx.*yy.*h2sr2;

e(:,:,3)=(-.25i*khi2)*  ( h0*(-s(3,:)/ep)+ xxh1sr*s(5,:)-yyh1sr*s(4,:)+(1i*gama/ep)*(xxh1sr*s(1,:)+yyh1sr*s(2,:)));
e(:,:,6)=(-.25i*khi2)*( h0*(-s(6,:)/mu)+xxh1sr*s(2,:)-yyh1sr*s(1,:)+(1i*gama/mu)*(xxh1sr*s(4,:)+yyh1sr*s(5,:)));
e(:,:,1)=(-.25i*mu)*( yyh1sr*((khi2/mu)*s(6,:))+h1sr*(-s(1,:)+(1i*gama/mu)*s(5,:))...
	-khi2*(xxyyh2sr2*s(2,:)-yy2h2sr2*s(1,:)+(1i*gama/mu)*(xxyyh2sr2*s(4,:)+yy2h2sr2*s(5,:))  ))...
	+(gama/4)* ( xxh1sr*((khi2/ep)*s(3,:))+h1sr*(s(5,:)+(1i*gama/ep)*s(1,:)) ...
	-khi2*(xx2h2sr2*s(5,:)-xxyyh2sr2*s(4,:)+(1i*gama/ep)*(xx2h2sr2*s(1,:)+xxyyh2sr2*s(2,:)))  );
e(:,:,2)=(gama/4)* ( yyh1sr*((khi2/ep)*s(3,:))+h1sr*(-s(4,:)+(1i*gama/ep)*s(2,:)) ...
	-khi2*(xxyyh2sr2*s(5,:)-yy2h2sr2*s(4,:)+(1i*gama/ep)*(xxyyh2sr2*s(1,:)+yy2h2sr2*s(2,:)))  )...
    +(.25i*mu)*( xxh1sr*((khi2/mu)*s(6,:))+h1sr*(s(2,:)+(1i*gama/mu)*s(4,:)) ...
	-khi2*(xx2h2sr2*s(2,:)-xxyyh2sr2*s(1,:)+(1i*gama/mu)*(xx2h2sr2*s(4,:)+xxyyh2sr2*s(5,:))  ));
e(:,:,4)=(-.25i*ep)*( yyh1sr*((khi2/ep)*s(3,:))+h1sr*(-s(4,:)+(1i*gama/ep)*s(2,:))...
	-khi2*(   xxyyh2sr2*s(5,:)-yy2h2sr2*s(4,:)+(1i*gama/ep)*(xxyyh2sr2*s(1,:)+yy2h2sr2*s(2,:))  ))...
	+(gama/4)* ( xxh1sr*((khi2/mu)*s(6,:))+h1sr*(s(2,:)+(1i*gama/mu)*s(4,:))...
	-khi2*(xx2h2sr2*s(2,:)-xxyyh2sr2*s(1,:)+(1i*gama/mu)*(xx2h2sr2*s(4,:)+xxyyh2sr2*s(5,:))  ));
e(:,:,5)=(gama/4)* ( yyh1sr*((khi2/mu)*s(6,:))+h1sr*(-s(1,:)+(1i*gama/mu)*s(5,:))...
	-khi2*(xxyyh2sr2*s(2,:)-yy2h2sr2*s(1,:)+(1i*gama/mu)*(xxyyh2sr2*s(4,:)+yy2h2sr2*s(5,:))  ))...
    +(.25i*ep)*( xxh1sr*((khi2/ep)*s(3,:))+h1sr*(s(5,:)+(1i*gama/ep)*s(1,:))...
	-khi2*(xx2h2sr2*s(5,:)-xxyyh2sr2*s(4,:)+(1i*gama/ep)*(xx2h2sr2*s(1,:)+xxyyh2sr2*s(2,:))  ));
e=permute(e,[1,3,2]);
if khi2<0;e(r==0,:,:)=0;end;% energie sur la source
if clone==0;e(:,4:6)=-1i*e(:,4:6);end;
if ~xmesh;e=reshape(e,length(z),length(x),length(y),6,size(s,2));end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function e=retpoint_012D(ep,mu,s,xs,x,y,z,clone) ;
% Il est possible de traiter simultanÈment plusieurs sources
% size(S)=[3,nb_Sources];en 1D   size(S)=[6,nb_Sources];en 1D
% alors e a une dimension de plus (nb_Sources)
if size(s,1)==1 | size(s,2)==1 ;s=s(:);end; 
%n=sqrt(ep*mu);
n=retsqrt(ep*mu,-1);% modif 2012
switch length(xs);
case 3;  % cas ' 2D ' (source dans espace 3 D )
if nargin<8;clone=0;end;xmesh=imag(clone)~=0;clone=real(clone);
if clone==1;s(1:3,:)=-i*s(1:3,:);end;
if xmesh;xx=x;yy=y;zz=z;else;[zz,xx,yy]=ndgrid(z,x,y);end;
xx=xx(:)-xs(1);yy=yy(:)-xs(2);zz=zz(:)-xs(3);


tri_PBA=0;if tri_PBA;% modif 12 2014 pour vitesse (probleme PBA)
[aaaa,kaaa,kkaaa]=retelimine([xx,yy,zz],1.e-10+i);%retcompare(aaaa(kkaaa,:),[xx,yy,zz])%
xx=aaaa(:,1);yy=aaaa(:,2);zz=aaaa(:,3);end;
r=sqrt(zz.^2+xx.^2+yy.^2);
%if isreal(xs);r=sqrt(zz.^2+xx.^2+yy.^2);else;r=-retsqrt(zz.^2+xx.^2+yy.^2);end;% pour faisceaux gaussiens
e=zeros(numel(r),size(s,2),6);
[eee,eeee,eeeee]=calfgh(n*r);eee=-(n/(4*pi))*eee;eeee=-(n^3/(4*pi))*eeee;eeeee=-(n^5/(4*pi))*eeeee;
%     ˘˘˘˘˘˘˘˘˘˘˘˘˘˘˘˘˘˘˘˘˘˘˘˘˘
%     e=zeros(numel(r),6);
%     [eee,eeee,eeeee]=calfgh(n*r);eee=-(n/(4*pi))*eee;eeee=-(n^3/(4*pi))*eeee;eeeee=-(n^5/(4*pi))*eeeee;
%     e(:,1)=-i*s(1)*(mu*eee+eeee/ep)-i*xx.*(s(1)*xx+s(2)*yy+s(3)*zz).*eeeee/ep-(yy*s(6)-zz*s(5)).*eeee;
%     e(:,2)=-i*s(2)*(mu*eee+eeee/ep)-i*yy.*(s(1)*xx+s(2)*yy+s(3)*zz).*eeeee/ep-(zz*s(4)-xx*s(6)).*eeee;
%     e(:,3)=-i*s(3)*(mu*eee+eeee/ep)-i*zz.*(s(1)*xx+s(2)*yy+s(3)*zz).*eeeee/ep-(xx*s(5)-yy*s(4)).*eeee;
%     e(:,4)=-(yy*s(3)-zz*s(2)).*eeee+i*s(4)*(ep*eee+eeee/mu)+i*xx.*(s(4)*xx+s(5)*yy+s(6)*zz).*eeeee/mu;
%     e(:,5)=-(zz*s(1)-xx*s(3)).*eeee+i*s(5)*(ep*eee+eeee/mu)+i*yy.*(s(4)*xx+s(5)*yy+s(6)*zz).*eeeee/mu;
%     e(:,6)=-(xx*s(2)-yy*s(1)).*eeee+i*s(6)*(ep*eee+eeee/mu)+i*zz.*(s(4)*xx+s(5)*yy+s(6)*zz).*eeeee/mu;
%     if clone==1;e(:,4:6)=i*e(:,4:6);end;
%     if ~xmesh e=reshape(e,length(z),length(x),length(y),6);end;
%     ˘˘˘˘˘˘˘˘˘˘˘˘˘˘˘˘˘˘˘˘˘˘˘˘˘˘˘˘˘˘˘
if (nnz(s)/numel(s)<.2 & tri_PBA);s=sparse(s);end;% pour gain de temps 12 2014
A1=1i*(mu*eee+eeee/ep);
Cx=xx.*eeee;Cy=yy.*eeee;Cz=zz.*eeee;
Cxx=xx.*xx.*eeeee;Cyy=yy.*yy.*eeeee;Czz=zz.*zz.*eeeee;Cxy=xx.*yy.*eeeee;Cxz=xx.*zz.*eeeee;Cyz=yy.*zz.*eeeee;
A2=1i*(ep*eee+eeee/mu);
e(:,:,1)=-A1*s(1,:)-(1i/ep)*(Cxx*s(1,:)+Cxy*s(2,:)+Cxz*s(3,:))-(Cy*s(6,:)-Cz*s(5,:));
e(:,:,2)=-A1*s(2,:)-(1i/ep)*(Cxy*s(1,:)+Cyy*s(2,:)+Cyz*s(3,:))-(Cz*s(4,:)-Cx*s(6,:));
e(:,:,3)=-A1*s(3,:)-(1i/ep)*(Cxz*s(1,:)+Cyz*s(2,:)+Czz*s(3,:))-(Cx*s(5,:)-Cy*s(4,:));
e(:,:,4)=-(Cy*s(3,:)-Cz*s(2,:))+A2*s(4,:)+(1i/mu)*(Cxx*s(4,:)+Cxy*s(5,:)+Cxz*s(6,:));
e(:,:,5)=-(Cz*s(1,:)-Cx*s(3,:))+A2*s(5,:)+(1i/mu)*(Cxy*s(4,:)+Cyy*s(5,:)+Cyz*s(6,:));
e(:,:,6)=-(Cx*s(2,:)-Cy*s(1,:))+A2*s(6,:)+(1i/mu)*(Cxz*s(4,:)+Cyz*s(5,:)+Czz*s(6,:));
e=permute(e,[1,3,2]);

if clone==1;e(:,4:6,:)=1i*e(:,4:6,:);end;
if tri_PBA;e=e(kkaaa,:,:);end;

if ~xmesh e=reshape(e,length(z),length(x),length(y),6,size(s,2));end;

case 2;  % cas ' 1D ' (source dans espace 2 D )
if nargin<8;clone=0;end;xmesh=imag(clone)~=0;clone=real(clone);
if clone==1;s(1,:)=-i*s(1,:);end;
if z==2;[ep,mu]=deal(mu,ep);end;% z=pol
if xmesh;xx=x;yy=y;else;[yy,xx]=ndgrid(y,x);end;
xx=xx(:)-xs(1);yy=yy(:)-xs(2);
r=sqrt(xx.^2+yy.^2);
%if isreal(xs);r=sqrt(xx.^2+yy.^2);else;r=-retsqrt(xx.^2+yy.^2);end;% pour faisceaux gaussiens
e=zeros(numel(r),size(s,2),3);
[h0,hy,hx,hyy,hxy,hxx]=calh(n*r,n*xx,n*yy);

e(:,:,1)=-(h0*s(1,:)*mu/4)+  (i*n/4)*(hy*s(2,:) -hx*s(3,:));
e(:,:,2)=-(hy*s(1,:)*n*i/4) -(n^2/(4*mu))*(hyy*s(2,:) -hxy*s(3,:));
e(:,:,3)=(hx*s(1,:)*n*i/4)+(n^2/(4*mu))*(hxy*s(2,:)-hxx*s(3,:));
e=permute(e,[1,3,2]);
if clone==1;e(:,2:3,:)=i*e(:,2:3,:);end;

if ~xmesh e=reshape(e,length(y),length(x),3,size(s,2));end;

case 1;  % cas ' 0D ' (source dans espace 1 D )
if nargin<7;clone=0;else;clone=real(z);end; 
if clone==1;s(1,:)=-i*s(1,:);end;
if y==2;muu=mu;mu=ep;ep=muu;end;
ap=(s(2,:)-s(1,:)*sqrt(mu/ep))/2;am=-(s(2,:)+s(1,:)*sqrt(mu/ep))/2;
x=x-xs;[f0,f]=retfind(x==0);[fp,fm]=retfind(x(f)>0);fp=f(fp);fm=f(fm);f0=reshape(f0,length(f0),1);fp=reshape(fp,length(fp),1);fm=reshape(fm,length(fm),1);% pour les vecteurs vides    
e=zeros(length(x),size(s,2),2);    
e(fp,:,1)=(exp(i*n*x(fp)))*ap;e(fm,:,1)=(exp(-i*n*x(fm)))*am;e(f0,:,1)=(ones(size(f0))/2)*(ap+am);
e(fp,:,2)=sqrt(ep/mu)*e(fp,:,1);e(fm,:,2)=-sqrt(ep/mu)*e(fm,:,1);e(f0,:,2)=(ones(size(f0))*sqrt(ep/mu)/2)*(ap-am);
e=permute(e,[1,3,2]);
if clone==1;e(:,2)=i*e(:,2);end;
end;    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [h0,hy,hx,hyy,hxy,hxx]=calh(r,x,y);
[h0,hy,hx,hyy,hxy,hxx]=deal(zeros(size(r)));
[g,f]=retfind(abs(r)<100*eps);
r=r(f);x=x(f);y=y(f);



if ~isempty(f);
% pour calculer correctement la partie reelle si r est petit
[ff,fff]=retfind(abs(imag(r))<.5);% modif 8 2013 pour r complexe
[h0r,h1r,h2r,h0i,h1i,h2i]=deal(zeros(size(r)));
h0r(ff)=retbessel('j',0,r(ff));h1r(ff)=retbessel('j',1,r(ff));h2r(ff)=retbessel('j',2,r(ff));
h0i(ff)=retbessel('y',0,r(ff));h1i(ff)=retbessel('y',1,r(ff));h2i(ff)=retbessel('y',2,r(ff));
h0r(fff)=retbessel('h',0,1,r(fff));h1r(fff)=retbessel('h',1,1,r(fff));h2r(fff)=retbessel('h',2,1,r(fff));

hyr=h1r.*y./r;hxr=h1r.*x./r;
hyyr=(y.^2./r.*h2r-h1r)./r;hxyr=x.*y./(r.^2).*h2r;hxxr=(x.^2./r.*h2r-h1r)./r;
hyi=h1i.*y./r;hxi=h1i.*x./r;
hyyi=(y.^2./r.*h2i-h1i)./r;hxyi=x.*y./(r.^2).*h2i;hxxi=(x.^2./r.*h2i-h1i)./r;
h0(f)=h0r+i*h0i;hy(f)=hyr+i*hyi;hx(f)=hxr+i*hxi;hyy(f)=hyyr+i*hyyi;hxy(f)=hxyr+i*hxyi;hxx(f)=hxxr+i*hxxi; 
end;

h0(g)=1;hy(g)=0;hx(g)=0;hyy(g)=-.5;hxy(g)=0;hxx(g)=-.5;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [h0,h1sr,h2sr2]=calh_conique(r);
[h0,h1sr,h2sr2]=deal(zeros(size(r)));[g,f]=retfind(abs(r)<100*eps);
[g,f]=retfind(abs(r)<100*eps);
r=r(f);
if ~isempty(f);
% pour calculer correctement la partie reelle si r est petit
[ff,fff]=retfind(abs(imag(r))<.5);% modif 8 2013 pour r complexe
[h0r,h1srr,h2sr2r,h0i,h1sri,h2sr2i]=deal(zeros(size(r)));
h0r(ff)=retbessel('j',0,r(ff));h1srr(ff)=retbessel('j',1,r(ff))./r(ff);h2sr2r(ff)=retbessel('j',2,r(ff))./(r(ff).^2);
h0i(ff)=retbessel('y',0,r(ff));h1sri(ff)=retbessel('y',1,r(ff))./r(ff);h2sr2i(ff)=retbessel('y',2,r(ff))./(r(ff).^2);
h0r(fff)=retbessel('h',0,1,r(fff));h1srr(fff)=retbessel('h',1,1,r(fff))./r(fff);h2sr2r(fff)=retbessel('h',2,1,r(fff))./(r(fff).^2);

h0(f)=h0r+1i*h0i;h1sr(f)=h1srr+1i*h1sri;h2sr2(f)=h2sr2r+1i*h2sr2i; 

end;
h0(g)=1;h1sr(g)=.5;h2sr2(g)=1/8; 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [eee,eeee,eeeee]=calfgh(r);
eee=zeros(size(r));eeee=zeros(size(r));eeeee=zeros(size(r));
[fff,g]=retfind(abs(r)<100*eps);

[f,ff]=retfind(abs(r(g))>1.e-2);f=g(f);ff=g(ff);
eee(f)=exp(i*r(f))./r(f);eeee(f)=(i-1./r(f)).*eee(f)./r(f);eeeee(f)=(-1-3i./r(f)+3./(r(f).^2)).*eee(f)./(r(f).^2);
%    -S0(r)                    -i/r S1(r)                                        i/r^2 S2(r) 
% Ce programme est 5 fois plus rapide que ret_bessel_spherique et aussi precis !

rr=r(ff).^2;
eee(ff)=polyval([1/40320,-1/720,1/24,-1/2,1],rr)./r(ff)+i*polyval([1/362880,-1/5040,1/120,-1/6,1],rr);
eeee(ff)=polyval([1/5760,-1/144,1/8,-1/2,-1],rr)./(r(ff).*r(ff).*r(ff))+i*polyval([-1/3991680,1/45360,-1/840,1/30,-1/3],rr);
eeeee(ff)=polyval([1/1152,-1/48,1/8,1/2,3],rr)./(r(ff).*(r(ff).*r(ff)).^2)+i*polyval([1/51891840,-1/498960,1/7560,-1/210,1/15],rr); 

% rr=r(fff).^2;
% eee(fff)=i*polyval([1/362880,-1/5040,1/120,-1/6,1],rr);
% eeee(fff)=i*polyval([-1/3991680,1/45360,-1/840,1/30,-1/3],rr);
% eeeee(fff)=i*polyval([1/51891840,-1/498960,1/7560,-1/210,1/15],rr);
eee(fff)=1i;
eeee(fff)=-1i/3;
eeeee(fff)=1i/15;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function e=retpoint_periodique_2D(ep,k0,S,centre,x,y,z,d,beta0,zmax,clone,parm);

Nrab=1;% diracs de chaque cote
if nargin<12;parm=[];end;if nargin<11;clone=0;end;if nargin<10;zmax=0;end;

%if isstruct(clone);[clone,parm]=deal(zmax,clone);zmax=0;end;% test 

if length(d)==1;e=retpoint_chaine_2D(ep,k0,S,centre,x,y,z,d,beta0,clone,parm);return;end;% chaine 

defaultopt=struct('point',1);
point=retoptimget(parm,'point',defaultopt,'fast');
n=retsqrt(ep/k0,-1);
%zmax=max([max(abs(z(:)-centre(3))),zmax,1.e-1*abs(real(2*pi/(k0*n)))]);
zmax=max([max(abs(z(:)-centre(3))),zmax,min(1.e-1*abs(real(2*pi/(k0*n))),min(d/5))]);% modif 1 2015
%zmax=max([max(abs(z(:)-centre(3))),zmax,1.e-2*abs(real(2*pi/(k0*n)))]);% modif 1 2 2014
init=retpoint_periodique_2D_store(d,n,beta0,k0,zmax,Nrab);

xmesh=imag(clone)~=0;clone=real(clone);
if ~xmesh;nx=length(x);ny=length(y);nz=length(z);[z,x,y]=ndgrid(z,x,y);end;
x=x(:)-centre(1);y=y(:)-centre(2);z=z(:)-centre(3);
p1=round(x/(2*d(1)));p2=round(y/(2*d(2)));xx=x-2*d(1)*p1;yy=y-2*d(2)*p2;% periodisation
PP=retelimine(p1+1i*p2);

if min(size(S))==1;S=S(:);end;nb_S=size(S,2);% S classe par colonnes
if clone==0;S(1:3,:)=i*S(1:3,:);end;% clonage de la source
e=retpoint_periodique_2D_interpol(init,xx,yy,z,S);

for ii=1:length(PP);f=find((p1+1i*p2)==PP(ii));if ~isempty(f);
if PP(ii)~=0;
offset_x=2*real(PP(ii))*d(1);offset_y=2*imag(PP(ii))*d(2);amp=exp(1i*beta0(1)*offset_x+1i*beta0(2)*offset_y);
for ix=-Nrab:Nrab;for iy=-Nrab:Nrab;% on rajoute les sources laterales et le centre dans les coordonnees shiftÈes (xx(f) ...
e(f,:,:)=e(f,:,:)+exp(1i*beta0(1)*ix*d(1)+1i*iy*beta0(2)*d(2))*retpoint(ep,k0,S,[ix*d(1),iy*d(2),0],xx(f),yy(f),z(f),1+i);      
end;end;
e(f,:,:)=amp*e(f,:,:)-retpoint(ep,k0,S,[0,0,0],x(f),y(f),z(f),1+i);% on enleve le centre dans les coordonnees non shiftÈes
else;
for ix=-Nrab:Nrab;for iy=-Nrab:Nrab;% on rajoute les sources laterales et par centre
if ix~=0 | iy~=0;e(f,:,:)=e(f,:,:)+exp(1i*beta0(1)*ix*d(1)+1i*iy*beta0(2)*d(2))*retpoint(ep,k0,S,[ix*d(1),iy*d(2),0],x(f),y(f),z(f),1+i);end;      
end;end;
end;end;end;% ii

%     for ix=-Nrab:Nrab;for iy=-Nrab:Nrab;% on rajoute les sources laterales
%     %if ix~=0 | iy~=0;for ii=1:nb_S;e(:,:,ii)=e(:,:,ii)+exp(1i*beta0(1)*ix*d(1)+1i*iy*beta0(2)*d(2))*retpoint(ep,k0,S(:,ii),[ix*d(1),iy*d(2),centre(3)],xx,yy,z,1+i);end;end;        
%     if ix~=0 | iy~=0;e=e+exp(1i*beta0(1)*ix*d(1)+1i*iy*beta0(2)*d(2))*retpoint(ep,k0,S,[ix*d(1),iy*d(2),0],x,y,z,1+i);end;       
%     end;end;

if point==1;e=e+retpoint(ep,k0,S,[0,0,0],x,y,z,1+i);end;% on ajoute la source non periodisÈe % modif 11 2012

% [offset_x,prv,ii]=retelimine(x(:)-centre(1)-xx(:));[offset_y,prv,jj]=retelimine(y(:)-centre(2)-yy(:));% modif 4 2013 pseudo periodisation
% if max(abs(offset_x))>100*eps*d(1) | max(abs(offset_y))>100*eps*d(2);amp_x=exp(1i*beta0(1)*offset_x);amp_y=exp(1i*beta0(2)*offset_y);for kk=1:nb_S;e(:,:,kk)=retdiag(amp_x(ii).*amp_y(jj))*e(:,:,kk);end;end;

if clone==0;e(:,4:6,:)=-1i*e(:,4:6,:);end;% declonage du champ
if ~xmesh;e=reshape(e,nz,nx,ny,6,nb_S);end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function e=retpoint_periodique_2D_sos(ep,k0,S,centre,x,y,z,d,beta0,zmax,clone,parm);

Nrab=1;% diracs de chaque cote
if nargin<12;parm=[];end;if nargin<11;clone=0;end;if nargin<10;zmax=0;end;

%if isstruct(clone);[clone,parm]=deal(zmax,clone);zmax=0;end;% test 

if length(d)==1;e=retpoint_chaine_2D(ep,k0,S,centre,x,y,z,d,beta0,clone,parm);return;end;% chaine 

defaultopt=struct('point',1);
point=retoptimget(parm,'point',defaultopt,'fast');
n=retsqrt(ep/k0,-1);
zmax=max([max(abs(z(:)-centre(3))),zmax,1.e-1*abs(real(2*pi/(k0*n)))]);
init=retpoint_periodique_2D_store(d,n,beta0,k0,zmax,Nrab);

xmesh=imag(clone)~=0;clone=real(clone);
if ~xmesh;nx=length(x);ny=length(y);nz=length(z);[z,x,y]=ndgrid(z,x,y);end;
xx=mod(x+d(1)-centre(1),2*d(1))-d(1);yy=mod(y+d(2)-centre(2),2*d(2))-d(2);% modif 4 2013

if min(size(S))==1;S=S(:);end;nb_S=size(S,2);% S classe par colonnes
if clone==0;S(1:3,:)=i*S(1:3,:);end;% clonage de la source
e=retpoint_periodique_2D_interpol(init,xx(:),yy(:),z(:)-centre(3),S);
for ix=-Nrab:Nrab;for iy=-Nrab:Nrab;% on rejoute les sources laterales
%if ix~=0 | iy~=0;for ii=1:nb_S;e(:,:,ii)=e(:,:,ii)+exp(1i*beta0(1)*ix*d(1)+1i*iy*beta0(2)*d(2))*retpoint(ep,k0,S(:,ii),[ix*d(1),iy*d(2),centre(3)],xx,yy,z,1+i);end;end;        
if ix~=0 | iy~=0;e=e+exp(1i*beta0(1)*ix*d(1)+1i*iy*beta0(2)*d(2))*retpoint(ep,k0,S,[ix*d(1),iy*d(2),centre(3)],xx,yy,z,1+i);end;       
end;end;

%if point==1;for ii=1:nb_S;e(:,:,ii)=e(:,:,ii)+retpoint(ep,k0,S(:,ii),[0,0,centre(3)],xx,yy,z,1+i);end;end;% on ajoute la source non periodisÈe % modif 11 2012
if point==1;e=e+retpoint(ep,k0,S,[0,0,centre(3)],xx,yy,z,1+i);end;% on ajoute la source non periodisÈe % modif 11 2012
[offset_x,prv,ii]=retelimine(x(:)-centre(1)-xx(:));[offset_y,prv,jj]=retelimine(y(:)-centre(2)-yy(:));% modif 4 2013 pseudo periodisation
if max(abs(offset_x))>100*eps*d(1) | max(abs(offset_y))>100*eps*d(2);amp_x=exp(1i*beta0(1)*offset_x);amp_y=exp(1i*beta0(2)*offset_y);for kk=1:nb_S;e(:,:,kk)=retdiag(amp_x(ii).*amp_y(jj))*e(:,:,kk);end;end;

if clone==0;e(:,4:6,:)=-1i*e(:,4:6,:);end;% declonage du champ
if ~xmesh;e=reshape(e,nz,nx,ny,6,nb_S);end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function e=retpoint_chaine_2D(ep,k0,S,centre,x,y,z,d,beta0,clone,parm);
defaultopt=struct('point',1);
point=retoptimget(parm,'point',defaultopt,'fast');
xmesh=imag(clone)~=0;clone=real(clone);
if ~xmesh;nx=length(x);ny=length(y);nz=length(z);[z,x,y]=ndgrid(z,x,y);end;
if min(size(S))==1;S=S(:);end;nb_S=size(S,2);% S classe par colonnes
if clone==0;S(1:3,:)=i*S(1:3,:);end;% clonage de la source
N=200;
[dd,XX]=ndgrid(-N:N,x);XXX=XX(:)-d*dd(:);
YY=repmat(y.',2*N+1,1);ZZ=repmat(z.',2*N+1,1);
e=zeros(length(XXX),6,nb_S);
for ii=1:nb_S;e(:,:,ii)=retpoint(ep,k0,S(:,ii),[0,0,0],XXX(:),YY(:),ZZ(:),1+i);end;
e=reshape(e,2*N+1,6*length(x)*nb_S);
if point~=1;e(N+1,:)=0;end; % pas de source au centre
apod=7.25;e=(retchamp([apod,2*N+1]).*exp(i*beta0(1)*(-N:N)*d))*e;
e=reshape(e,length(z),6,nb_S);
if clone==0;e(:,4:6,:)=-1i*e(:,4:6,:);end;% declonage du champ
if ~xmesh;e=reshape(e,nz,nx,ny,6,nb_S);end;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function varargout=retpoint_periodique_2D_store(varargin);% pour eviter de recommencer
persistent store_in store_out
if isempty(store_in);store_in={};store_out={};end;
%for ii=1:length(store_in);if isequal(varargin(1:end-1),store_in{ii}(1:end-1)) & varargin{end}<=store_in{ii}{end} & (length(store_out{ii})==nargout);varargout=store_out{ii};return;end;end;
for ii=1:length(store_in);if isequal(varargin(1:end-1),store_in{ii}(1:end-1)) & varargin{end}<=store_in{ii}{end} ;varargout=store_out{ii};return;end;end;% modif 12 2012
[varargout{1:nargout}]=retpoint_periodique_2D_init(varargin{:});
store_in=[store_in,{varargin}];store_out=[store_out,{varargout}];% on stocke le resultat

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function init=retpoint_periodique_2D_init(d,n,beta0,k0,zmax,Nrab);
pqr=[[0,0,0] ;[1,0,0] ;[0,1,0] ;[0,0,1] ;[2,0,0] ;[0,2,0] ;[0,0,2] ;[1,1,0] ;[1,0,1] ;[0,1,1]];

k0n=k0*n;
nzmin=7;
nbz=max(nzmin,ceil(9*zmax/abs(real(2*pi/k0n))));% minimum 5
z=linspace(0,zmax,nbz+1);z=z(2:end);% obligatoirement regulierement espacÈs
NN_min=64;
NNx=NN_min;
while NNx<=2^11;
u=beta0(1)+2*pi/d(1)*[-NNx,NNx-1];
ermax=max(abs(u).^2.*exp(-(z(1)*abs(u))));   
if ermax<1.e-20;break;end;
NNx=2*NNx;
end;
NNy=NN_min;
while NNy<=2^11
v=beta0(2)+2*pi/d(2)*[-NNy,NNy-1];
ermax=max(abs(v).^2.*exp(-(z(1)*abs(v))));   
if ermax<1.e-20;break;end;
NNy=2*NNy;
end;

nnx=ceil(NNx/NN_min);nny=ceil(NNy/NN_min);
u=beta0(1)+2*pi/d(1)*(-NNx:NNx-1);x=(d(1)/(2*NNx))*(-2*NNx:nnx:2*NNx);
v=beta0(2)+2*pi/d(2)*(-NNy:NNy-1);y=(d(2)/(2*NNy))*(-2*NNy:nny:2*NNy);
[U,V]=ndgrid(u,v);
khi=retbidouille(retsqrt(k0n^2-(U.^2+V.^2),-1));
K=(exp(1i*beta0(1)*x.')*exp(1i*beta0(2)*y));


[X,Y,R2,numR,GU,GV,GW]=deal(cell(2*Nrab+1,2*Nrab+1));
for ii=-Nrab:Nrab;for jj=-Nrab:Nrab;
[X{ii+Nrab+1,jj+Nrab+1},Y{ii+Nrab+1,jj+Nrab+1}]=ndgrid(x-ii*d(1),y-jj*d(2));[R2{ii+Nrab+1,jj+Nrab+1},prv,numR{ii+Nrab+1,jj+Nrab+1}]=retelimine(X{ii+Nrab+1,jj+Nrab+1}.^2+Y{ii+Nrab+1,jj+Nrab+1}.^2,10*eps);% pour gain de temps
end;end;
GGG=zeros(length(x),length(y),length(z),10);
for iz=1:length(z);Z=z(iz);

while nnx>1;
ermax=max([abs(u(NNx/2+1))^2*exp(-(Z*abs(u(NNx/2+1)))),abs(u(NNx/2+NNx))^2*exp(-(Z*abs(u(NNx/2+NNx))))]);   
if ermax>1.e-20;break;end;
nnx=nnx/2;NNx=NNx/2;u=u(NNx+1:NNx+2*NNx);khi=khi(NNx+1:NNx+2*NNx,:);
end
while nny>1;
ermax=max([abs(v(NNy/2+1))^2*exp(-(Z*abs(v(NNy/2+1)))),abs(v(NNy/2+NNy))^2*exp(-(Z*abs(v(NNy/2+NNy))))]);   
if ermax>1.e-20;break;end;
nny=nny/2;NNy=NNy/2;v=v(NNy+1:NNy+2*NNy);khi=khi(:,NNy+1:NNy+2*NNy);
end   

for ii=-Nrab:Nrab;for jj=-Nrab:Nrab;
RR=sqrt(R2{ii+Nrab+1,jj+Nrab+1}+Z^2);
GU{ii+Nrab+1,jj+Nrab+1}=exp(1i*k0n*RR)./(-4*pi*RR);
GV{ii+Nrab+1,jj+Nrab+1}=(1i*k0n./RR-1./(RR.^2)).*GU{ii+Nrab+1,jj+Nrab+1};
GW{ii+Nrab+1,jj+Nrab+1}=(-k0n^2./(RR.^2)-3i*k0n./(RR.^3)+3./(RR.^4)).*GU{ii+Nrab+1,jj+Nrab+1};
GU{ii+Nrab+1,jj+Nrab+1}=reshape(GU{ii+Nrab+1,jj+Nrab+1}(numR{ii+Nrab+1,jj+Nrab+1}),size(X{ii+Nrab+1,jj+Nrab+1}));GV{ii+Nrab+1,jj+Nrab+1}=reshape(GV{ii+Nrab+1,jj+Nrab+1}(numR{ii+Nrab+1,jj+Nrab+1}),size(X{ii+Nrab+1,jj+Nrab+1}));GW{ii+Nrab+1,jj+Nrab+1}=reshape(GW{ii+Nrab+1,jj+Nrab+1}(numR{ii+Nrab+1,jj+Nrab+1}),size(X{ii+Nrab+1,jj+Nrab+1}));
end;end;


FF=(1/(2i*prod(d)))*exp(1i*khi*abs(Z));
FFF=FF./khi;
% disp(rettexte(NNx,NNy,nnx,nny));
for ipqr=1:size(pqr,1);
switch ipqr;
case 1;F=FFF;	
case 2;F=1i*(retdiag(u)*FFF);
case 3;F=1i*(FFF*retdiag(v));
case 4;F=1i*FF;
case 5;F=-(retdiag(u.^2)*FFF);
case 6;F=-(FFF*retdiag(v.^2));
case 7;F=-FF.*khi;
case 8;F=-retdiag(u)*FFF*retdiag(v);
case 9;F=-(retdiag(u)*FF);
case 10;F=-(FF*retdiag(v));
end;
F=(4*NNx*NNy)*(ifft2(ifftshift(F)));
GGG(:,:,iz,ipqr)=K.*F([1:nnx:2*NNx,1:nnx:2*NNx,1],[1:nny:2*NNy,1:nny:2*NNy,1]);clear F;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for ii=-Nrab:Nrab;for jj=-Nrab:Nrab;A=exp(1i*beta0(1)*ii*d(1)+1i*beta0(2)*jj*d(2));
switch ipqr;
case 1;GGG(:,:,iz,ipqr)=GGG(:,:,iz,ipqr)-A*GU{ii+Nrab+1,jj+Nrab+1};	
case 2;GGG(:,:,iz,ipqr)=GGG(:,:,iz,ipqr)-A*X{ii+Nrab+1,jj+Nrab+1}.*GV{ii+Nrab+1,jj+Nrab+1};
case 3;GGG(:,:,iz,ipqr)=GGG(:,:,iz,ipqr)-A*Y{ii+Nrab+1,jj+Nrab+1}.*GV{ii+Nrab+1,jj+Nrab+1};
case 4;GGG(:,:,iz,ipqr)=GGG(:,:,iz,ipqr)-A*Z.*GV{ii+Nrab+1,jj+Nrab+1};
case 5;GGG(:,:,iz,ipqr)=GGG(:,:,iz,ipqr)-A*GV{ii+Nrab+1,jj+Nrab+1}-A*X{ii+Nrab+1,jj+Nrab+1}.^2.*GW{ii+Nrab+1,jj+Nrab+1};
case 6;GGG(:,:,iz,ipqr)=GGG(:,:,iz,ipqr)-A*GV{ii+Nrab+1,jj+Nrab+1}-A*Y{ii+Nrab+1,jj+Nrab+1}.^2.*GW{ii+Nrab+1,jj+Nrab+1};
case 7;GGG(:,:,iz,ipqr)=GGG(:,:,iz,ipqr)-A*GV{ii+Nrab+1,jj+Nrab+1}-A*Z.^2.*GW{ii+Nrab+1,jj+Nrab+1};
case 8;GGG(:,:,iz,ipqr)=GGG(:,:,iz,ipqr)-A*X{ii+Nrab+1,jj+Nrab+1}.*Y{ii+Nrab+1,jj+Nrab+1}.*GW{ii+Nrab+1,jj+Nrab+1};
case 9;GGG(:,:,iz,ipqr)=GGG(:,:,iz,ipqr)-A*X{ii+Nrab+1,jj+Nrab+1}.*Z.*GW{ii+Nrab+1,jj+Nrab+1};
case 10;GGG(:,:,iz,ipqr)=GGG(:,:,iz,ipqr)-A*Y{ii+Nrab+1,jj+Nrab+1}.*Z.*GW{ii+Nrab+1,jj+Nrab+1};
end;
end;end;% ii jj
end;end;% ipqr iz
clear F FF FFF GU GV GW U V khi;

%%% extrapolation au centre
G0=zeros(length(x),length(y),1,10);
% a1=z(2)^2/(z(2)^2-z(1)^2);a2=-z(1)^2/(z(2)^2-z(1)^2);
% G0(:,:,:,[1,2,3,5,6,7,8])=a1*GGG(:,:,1,[1,2,3,5,6,7,8])+a2*GGG(:,:,2,[1,2,3,5,6,7,8]);% les autres composantes sont nulles

%deg=min(5,length(z));Mat=repmat(z(:).^2,1,deg).^repmat(0:deg-1,length(z),1)\eye(length(z));
%Mat=zeros(length(z));for jj=1:length(z);Mat(:,jj)=cos(2*(jj-1)*acos(z.'/z(end)));end;Mat=(-1).^[0:length(z)-1]/Mat;% polynomes de Tchebitchev
%MMat=zeros(length(z));for jj=1:length(z);MMat(:,jj)=cos(2*(jj-1)*acos(z.'/z(end)));end;MMat=inv(MMat);Mat=(-1).^[0:length(z)-1]*MMat;% polynomes de Tchebitchev
%prvMat=Mat;prvMMat=MMat;



[MMat,Mat]=retMathieu(nbz);
% if nbz<=50;load retMathieu;MMat=MMat{nbz};
% else;
% MMat=zeros(nbz);for jj=1:nbz;MMat(:,jj)=cos(2*(jj-1)*acos(z(1:nbz).'/z(nbz)));end;MMat=inv(MMat);% polynomes de Tchebitchev
% end;
% Mat=(-1).^[0:nbz-1]*MMat;





%for ipqr=[1,2,3,5,6,7,8];for ii=1:length(z);G0(:,:,1,ipqr)=G0(:,:,1,ipqr)+GGG(:,:,ii,ipqr)*Mat(1,ii);end;end;

GGG=reshape(GGG,length(x)*length(y),length(z),10);
for ipqr=[4,9,10];GGG(:,:,ipqr)=GGG(:,:,ipqr)*retdiag(1./z);end;% on divise par z les impairs pour les rendre pairs
GGG=reshape(GGG,length(x),length(y),length(z),10);
%prvG0=G0;for ipqr=1:10;for ii=1:length(z);prvG0(:,:,1,ipqr)=prvG0(:,:,1,ipqr)+GGG(:,:,ii,ipqr)*prvMat(1,ii);end;end;
for ipqr=1:10;for ii=1:nbz;G0(:,:,1,ipqr)=G0(:,:,1,ipqr)+GGG(:,:,ii,ipqr)*Mat(1,ii);end;end;

% p=20;t=linspace(-1,1,2*p);t=t(p+1:2*p);Tch=zeros(p);for jj=1:p;Tch(:,jj)=cos(2*(jj-1)*acos(t.'));end;M=Tch\eye(p);
% symetrisation en z
parite=mod(pqr(:,end),2);

% for ii=1:10;fx=find(abs(x)<max(x)/4);fy=find(abs(y)<max(y)/4);retfig;figure;retcolor(y(fy),x(fx),real(GGG(fx,fy,1,ii)),2i);hold on;;pause;end;stop


N=2;Ix0=retminabs(x);Iy0=retminabs(y);fx=Ix0-N:Ix0+N;fy=Iy0-N:Iy0+N;
[G00,G10,G01,G20,G11,G02]=deal(zeros(nbz,length(parite)));
for ipqr=1:length(parite);for iz=1:nbz; 
GGGx=retderivee(x(fx),GGG(fx,fy,iz,ipqr),1);GGGy=retderivee(y(fy),GGG(fx,fy,iz,ipqr),2);
GGGxx=retderivee(x(fx),GGGx,1);GGGxy=.5*(retderivee(x(fx),GGGy,1)+retderivee(y(fy),GGGx,2));GGGyy=retderivee(y(fy),GGGy,2);
G00(iz,ipqr)=GGG(Ix0,Iy0,iz,ipqr);G10(iz,ipqr)=GGGx(N+1,N+1);G01(iz,ipqr)=GGGy(N+1,N+1);G11(iz,ipqr)=GGGxy(N+1,N+1);G20(iz,ipqr)=GGGxx(N+1,N+1);G02(iz,ipqr)=GGGyy(N+1,N+1);
end;end;
z=[0,z];
GGG=cat(3,G0,GGG);

init={k0,n,x,y,z,GGG,parite,MMat,G00,G10,G01,G20,G11,G02};


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function e=retpoint_periodique_2D_interpol(init,xi,yi,zi,S);
xi=xi(:);yi=yi(:);zi=zi(:);
[k0,n,x,y,z,GGG,parite,MMat,G00,G10,G01,G20,G11,G02]=deal(init{:});nbz=size(MMat,1);
GGGG=zeros(length(xi),length(parite));

dx=3*mean(diff(x));dy=3*mean(diff(y));dz=3*mean(diff(z));

[zzi,kk,kkk]=retelimine(zi);
for ii=1:length(zzi);%zzi(ii)
mmat=(cos(2*(0:nbz-1)*acos(zzi(ii)/z(nbz+1)))*MMat);
if abs(zzi(ii))<z(2);g00=mmat*G00;g10=mmat*G10;g01=mmat*G01;g20=mmat*G20;g11=mmat*G11;g02=mmat*G02;end;

    ffii=find(kkk==ii);
    ffx=retelimine(round((xi(ffii)-x(1))/(x(2)-x(1))))+1;ffx=retelimine([ffx+1,ffx-1,ffx+2,ffx-2,ffx+3,ffx-3,ffx]);ffx=ffx(ffx>0 & ffx<=length(x));
    ffy=retelimine(round((yi(ffii)-y(1))/(y(2)-y(1))))+1;ffy=retelimine([ffy+1,ffy-1,ffy+2,ffy-2,ffy+3,ffy-3,ffy]);ffy=ffy(ffy>0 & ffy<=length(y));
    fxb=[0;find(diff(ffx)>1);length(ffx)];fyb=[0;find(diff(ffy)>1);length(ffy)];% decomposition en blocs
    for iix=1:length(fxb)-1;fx=ffx(fxb(iix)+1:fxb(iix+1));for iiy=1:length(fyb)-1;fy=ffy(fyb(iiy)+1:fyb(iiy+1));
    fii= ffii( (xi(ffii)>=x(fx(1))) &  (xi(ffii)<=x(fx(end))) &  (yi(ffii)>=y(fy(1))) & (yi(ffii)<=y(fy(end))));
    if ~isempty(fii);    
        
Nx=length(fx);Ny=length(fy);
G=reshape(permute(GGG(fx,fy,2:nbz+1,:),[1,2,4,3]),[],nbz)*(mmat).';
G=reshape(G,length(fx),length(fy),[]);
%for ipqr=1:length(parite);GGGG(:,ipqr)=interpn(x(fx),y(fy),G(:,:,ipqr),xi,yi,'spleen');if parite(ipqr)==1;GGGG(:,ipqr)=GGGG(:,ipqr).*zzi;end;end;% pour gain de temps 4 2013
%for ipqr=1:length(parite);GGGG(:,ipqr)=interpn(x(fx),y(fy),real(G(:,:,ipqr)),xi,yi,'pchip')+1i*interpn(x(fx),y(fy),imag(G(:,:,ipqr)),xi,yi,'pchip');if parite(ipqr)==1;GGGG(:,ipqr)=GGGG(:,ipqr).*zzi;end;end;
for ipqr=1:length(parite);
    if abs(zzi(ii))<z(2);% pour les tres pres on utilise le developpement ‡ l'origine
gg=g00(ipqr)+g10(ipqr)*xi(fii)+g01(ipqr)*yi(fii)+.5*g20(ipqr)*xi(fii).^2+g11(ipqr)*xi(fii).*yi(fii)+.5*g02(ipqr)*yi(fii).^2;
g=g00(ipqr)*ones(Nx,Ny)+g10(ipqr)*x(fx).'*ones(1,Ny)+g01(ipqr)*ones(Nx,1)*y(fy)+.5*g20(ipqr)*x(fx).^2.'*ones(1,Ny)+g11(ipqr)*x(fx).'*y(fy)+.5*g02(ipqr)*ones(Nx,1)*y(fy).^2;
    else;gg=0;g=0;end;
    if length(fx)<20 & length(fy)<20;
  GGGG(fii,ipqr)=gg+ret_interp2(x(fx),y(fy),G(:,:,ipqr)-g,xi(fii),yi(fii));
    else;
  GGGG(fii,ipqr)=gg+interpn(x(fx),y(fy),G(:,:,ipqr)-g,xi(fii),yi(fii),'spleen');
    end
if parite(ipqr)==1;GGGG(fii,ipqr)=GGGG(fii,ipqr).*zzi(ii);end;
end;

%for kk=1:10;retfig;retcolor(y(fy),x(fx),real(G(:,:,kk)*(zzi(ii)^parite(kk))),2i);hold on;plot3(yi(fii),xi(fii),real(GGGG(fii,kk)),'*k');pause;end;%stop
    end;end;end;

end;

e=zeros(length(xi),6,size(S,2));
mu=k0;ep=k0*n^2;
for i_S=1:size(S,2);
e(:,1,i_S)=-mu*S(1,i_S)*GGGG(:,1)-(1/ep)*(S(1,i_S)*GGGG(:,5)+S(2,i_S)*GGGG(:,8)+S(3,i_S)*GGGG(:,9))-S(6,i_S)*GGGG(:,3)+S(5,i_S)*GGGG(:,4);    
e(:,2,i_S)=-mu*S(2,i_S)*GGGG(:,1)-(1/ep)*(S(1,i_S)*GGGG(:,8)+S(2,i_S)*GGGG(:,6)+S(3,i_S)*GGGG(:,10))-S(4,i_S)*GGGG(:,4)+S(6,i_S)*GGGG(:,2);        
e(:,3,i_S)=-mu*S(3,i_S)*GGGG(:,1)-(1/ep)*(S(1,i_S)*GGGG(:,9)+S(2,i_S)*GGGG(:,10)+S(3,i_S)*GGGG(:,7))-S(5,i_S)*GGGG(:,2)+S(4,i_S)*GGGG(:,3);     

e(:,4,i_S)=-ep*S(4,i_S)*GGGG(:,1)-(1/mu)*(S(4,i_S)*GGGG(:,5)+S(5,i_S)*GGGG(:,8)+S(6,i_S)*GGGG(:,9))-S(3,i_S)*GGGG(:,3)+S(2,i_S)*GGGG(:,4);         
e(:,5,i_S)=-ep*S(5,i_S)*GGGG(:,1)-(1/mu)*(S(4,i_S)*GGGG(:,8)+S(5,i_S)*GGGG(:,6)+S(6,i_S)*GGGG(:,10))-S(1,i_S)*GGGG(:,4)+S(3,i_S)*GGGG(:,2);    
e(:,6,i_S)=-ep*S(6,i_S)*GGGG(:,1)-(1/mu)*(S(4,i_S)*GGGG(:,9)+S(5,i_S)*GGGG(:,10)+S(6,i_S)*GGGG(:,7))-S(2,i_S)*GGGG(:,2)+S(1,i_S)*GGGG(:,3);     
end;
%figure;retcolor(yy,xx,real(ZZZ),2i);hold on;plot3(Y,X,real(Z),'*k');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function ZZZ=ret_interp2(x,y,Z,XX,YY);

persistent aa;if isempty(aa); aa=cell(20,20);end;

Nx=length(x);Ny=length(y);
%x0=mean(x);y0=mean(y);xn=x-x0;yn=y-y0;Ix=max(abs(xn));Iy=max(abs(yn));
x0=(x(1)+x(end))/2;y0=(y(1)+y(end))/2;xn=x-x0;yn=y-y0;Ix=abs(xn(1));Iy=abs(yn(1));
if Nx>size(aa,1) | Ny>size(aa,2) |  isempty(aa{Nx,Ny});
xn=xn/Ix;yn=yn/Iy;
Xnp=cos( repmat(acos(xn(:)),1,Nx).*repmat([0:Nx-1],numel(xn),1));
Ynq=cos( repmat(acos(yn(:)),1,Ny).*repmat([0:Ny-1],numel(yn),1));
%a=reshape( (reshape(permute(reshape(Xnp(:)*Ynq(:).',length(xn),Nx,length(yn),Ny),[1,3,2,4]),Nx*Ny,Nx*Ny))\Z(:),Nx,Ny);
a=reshape(permute(reshape(Xnp(:)*Ynq(:).',length(xn),Nx,length(yn),Ny),[1,3,2,4]),Nx*Ny,Nx*Ny);
a=inv(a);
aa{Nx,Ny}=a;
else;
a=aa{Nx,Ny};
end

Z0=mean(Z(:));
a=reshape(a*(Z(:)-Z0),Nx,Ny);
B=zeros(numel(XX),Nx,Ny);
XXn=(XX-x0)/Ix;
YYn=(YY-y0)/Iy;

ZZZ=0;
XXnp=cos( repmat(acos(XXn(:)),1,Nx).*repmat([0:Nx-1],numel(XXn),1));
YYnq=cos( repmat(acos(YYn(:)),1,Ny).*repmat([0:Ny-1],numel(YYn),1));

for iq=0:Ny-1;
for ip=0:Nx-1;
ZZZ=ZZZ+a(ip+1,iq+1)*(XXnp(:,ip+1).*YYnq(:,iq+1));
end;end;
ZZZ=Z0+reshape(ZZZ,size(XX));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function e=retpoint_periodique_1D(ep,k0,S,centre,x,y,pol,d,beta0,ymax,clone,parm);
Nrab=1;% diracs de chaque cote
if nargin<12;parm=[];end;if nargin<11;clone=0;end;if nargin<10;ymax=0;end;
defaultopt=struct('point',1,'Ewald',0);
point=retoptimget(parm,'point',defaultopt,'fast');
Ewald=retoptimget(parm,'Ewald',defaultopt,'fast');
Ewald=Ewald & numel(x)<100;
n=retsqrt(ep/k0,-1);
if  ~Ewald;
ymax=max([max(abs(y(:)-centre(2))),ymax,1.e-1*abs(real(2*pi/(k0*n)))]);
init=retpoint_periodique_1D_store(d,n,beta0,k0,ymax,Nrab);
end;
xmesh=imag(clone)~=0;clone=real(clone);
if ~xmesh;nx=length(x);ny=length(y);[y,x]=ndgrid(y,x);end;
x=x(:)-centre(1);y=y(:)-centre(2);
p=round(x/(2*d));xx=x-2*d*p;% periodisation
PP=retelimine(p);


if min(size(S))==1;S=S(:);end;nb_S=size(S,2);% S classe par colonnes
if clone==0;S(1,:)=1i*S(1,:);end;% clonage de la source
if Ewald;
e=ret_Formule_Petit(xx,y,S,d,n,pol,beta0,k0,Nrab+1);
else;
e=retpoint_periodique_1D_interpol(init,xx,y,pol,S);
end;

for ii=1:length(PP);f=find(p==PP(ii));if ~isempty(f);
if PP(ii)~=0;
offset=2*PP(ii)*d(1);amp=exp(1i*beta0*offset);
for ix=-Nrab:Nrab;% on rajoute les sources laterales et le centre dans les coordonnees shiftÈes (xx(f) ...
e(f,:,:)=e(f,:,:)+exp(1i*beta0*ix*d)*retpoint(ep,k0,S,[ix*d,0],xx(f),y(f),pol,1+1i);      
end;
e(f,:,:)=amp*e(f,:,:)-retpoint(ep,k0,S,[0,0],x(f),y(f),pol,1+1i);
else;
for ix=-Nrab:Nrab;% on rajoute les sources laterales et pas centre
if ix~=0 ;e(f,:,:)=e(f,:,:)+exp(1i*beta0*ix*d)*retpoint(ep,k0,S,[ix*d,0],x(f),y(f),pol,1+1i);end;      
end;
end;end;end;% ii

%     for ix=-Nrab:Nrab;for iy=-Nrab:Nrab;% on rajoute les sources laterales
%     %if ix~=0 | iy~=0;for ii=1:nb_S;e(:,:,ii)=e(:,:,ii)+exp(1i*beta0(1)*ix*d(1)+1i*iy*beta0(2)*d(2))*retpoint(ep,k0,S(:,ii),[ix*d(1),iy*d(2),centre(3)],xx,yy,z,1+i);end;end;        
%     if ix~=0 | iy~=0;e=e+exp(1i*beta0(1)*ix*d(1)+1i*iy*beta0(2)*d(2))*retpoint(ep,k0,S,[ix*d(1),iy*d(2),0],x,y,z,1+i);end;       
%     end;end;

%if point==1;for ii=1:nb_S;e(:,:,ii)=e(:,:,ii)+retpoint(ep,k0,S(:,ii),[0,centre(2)],xx,y,pol,1+i);end;end;% on ajoute la source non periodisÈe % modif 11 2012
%if point==1;e=e+retpoint(ep,k0,S,[0,centre(2)],xx,y,pol,1+1i);end;% on ajoute la source non periodisÈe % modif 11 2012
if point==1;e=e+retpoint(ep,k0,S,[0,0],x,y,pol,1+1i);end;% on ajoute la source non periodisÈe % modif 11 2012
% [offset,prv,ii]=retelimine(x(:)-xx(:));% modif 4 2013 pseudo periodisation
% if max(abs(offset))>100*eps*d ;amp=exp(1i*beta0*offset);e=retdiag(amp(ii))*e;end;

if clone==0;e(:,2:end,:)=-1i*e(:,2:end,:);end;% declonage du champ
if ~xmesh;e=reshape(e,ny,nx,[],nb_S);end;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function varargout=retpoint_periodique_1D_store(varargin);% pour eviter de recommencer
persistent store_in store_out
if isempty(store_in);store_in={};store_out={};end;
for ii=1:length(store_in);if isequal(varargin(1:end-1),store_in{ii}(1:end-1)) & varargin{end}<=store_in{ii}{end} & (length(store_out{ii})==nargout);varargout=store_out{ii};return;end;end;
[varargout{1:nargout}]=retpoint_periodique_1D_init(varargin{:});
store_in=[store_in,{varargin}];store_out=[store_out,{varargout}];% on stocke le resultat

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function init=retpoint_periodique_1D_init(d,n,beta0,k0,ymax,Nrab);
pq=[[0,0] ;[1,0] ;[0,1] ;[2,0,] ;[0,2]  ;[1,1] ];

k0n=k0*n;nymin=7;
nby=max(nymin,ceil(9*ymax/abs(real(2*pi/k0n))));%minimum 5
y=linspace(0,ymax,nby+1);y=y(2:end);% obligatoirement regulierement espacÈs
NN_min=128;
NN=NN_min;
while NN<=2^12
u=beta0+2*pi/d*[-NN,NN-1];
ermax=max(abs(u).^2.*exp(-(y(1)*abs(u))));   
if ermax<1.e-20;break;end;
NN=2*NN;
end;
nn=ceil(NN/NN_min);
u=beta0+(2*pi/d)*(-NN:NN-1);x=(d/(2*NN))*(-2*NN:nn:2*NN);u=u(:);x=x(:);
khi=retbidouille(retsqrt(k0n^2-u.^2,-1));
K=exp(1i*beta0*x);
GGG=zeros(length(x),length(y),size(pq,1));
for iy=1:length(y);% iy ‡‡‡‡‡‡‡‡‡‡‡‡‡‡‡‡‡‡‡‡‡‡‡
while nn>1;
ermax=max([abs(u(NN/2+1))^2*exp(-(y(1)*abs(u(NN/2+1)))),abs(u(NN/2+NN))^2*exp(-((1)*abs(u(NN/2+NN))))]);   
if ermax>1.e-20;break;end;
nn=nn/2;NN=NN/2;u=u(NN+1:NN+2*NN);khi=khi(NN+1:NN+2*NN);
end


for ii=-Nrab:Nrab;
R=sqrt((x-ii*d).^2+y(iy).^2);
GU{ii+Nrab+1}=retbessel('h',[0,1,2],1,k0n*R);
GV{ii+Nrab+1}=-k0n*GU{ii+Nrab+1}(:,2)./(4i*R);
GW{ii+Nrab+1}=k0n^2*GU{ii+Nrab+1}(:,3)./(4i*R.^2);
GU{ii+Nrab+1}=GU{ii+Nrab+1}(:,1)/4i;
end;

FF=(1/(2i*d))*exp(1i*khi*abs(y(iy)));
FFF=FF./khi;
for ipq=1:size(pq,1);
switch ipq;
case 1;F=FFF;	
case 2;F=1i*(u.*FFF);
case 3;F=1i*FF;
case 4;F=-(u.^2.*FFF);
case 5;F=-FF.*khi;
case 6;F=-(u.*FF);
end;
F=(2*NN)*(ifft(ifftshift(F)));
GGG(:,iy,ipq)=K.*F([1:nn:2*NN,1:nn:2*NN,1]);clear F;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for ii=-Nrab:Nrab;A=exp(1i*beta0*ii*d);
switch ipq;
case 1;GGG(:,iy,ipq)=GGG(:,iy,ipq)-A*GU{ii+Nrab+1};	
case 2;GGG(:,iy,ipq)=GGG(:,iy,ipq)-A*(x-ii*d).*GV{ii+Nrab+1};
case 3;GGG(:,iy,ipq)=GGG(:,iy,ipq)-A*y(iy).*GV{ii+Nrab+1};
case 4;GGG(:,iy,ipq)=GGG(:,iy,ipq)-A*GV{ii+Nrab+1}-A*((x-ii*d).^2).*GW{ii+Nrab+1};
case 5;GGG(:,iy,ipq)=GGG(:,iy,ipq)-A*GV{ii+Nrab+1}-A*(y(iy).^2).*GW{ii+Nrab+1};
case 6;GGG(:,iy,ipq)=GGG(:,iy,ipq)-A*(x-ii*d).*y(iy).*GW{ii+Nrab+1};
end;
end; % ii
end;% ipq 
end;% iy ‡‡‡‡‡‡‡‡‡‡‡‡‡‡‡‡‡‡‡‡‡‡‡‡‡

clear F FF FFF GU GV GW U V khi;
%%% extrapolation au centre
G0=zeros(length(x),1,6);
[MMat,Mat]=retMathieu(nby);

% if nby<=50;load retMathieu;MMat=MMat{nby};
% else;
% MMat=zeros(nby);for jj=1:nby;MMat(:,jj)=cos(2*(jj-1)*acos(y(1:nby).'/y(nby)));end;MMat=inv(MMat);% polynomes de Tchebitchev
% end;
% Mat=(-1).^[0:nby-1]*MMat;

for ipq=[3,6];GGG(:,:,ipq)=GGG(:,:,ipq)*retdiag(1./y);end;% on divise par y les impairs pour les rendre pairs
for ipq=1:6;G0(:,1,ipq)=GGG(:,:,ipq)*Mat(1,:).';end;
% [prv,GGGG]=ret_Formule_Petit(x,0*ones(size(x)),[],d,n,0,beta0,k0,Nrab+1);G0=reshape(GGGG(:,[1,2,5,4,5,7]),length(x),1,6);% 10 2013
y=[0,y];
GGG=cat(2,G0,GGG);
parite=mod(pq(:,end),2);
init={k0,n,x,y,GGG,parite,MMat};


% for ii=1:6;fx=find(abs(x)<max(x)/4);retfig;figure;retcolor(y,x(fx),real(GGG(fx,:,ii)),2i);hold on;;pause;end;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function e=retpoint_periodique_1D_interpol(init,xi,yi,pol,S);
xi=xi(:);yi=yi(:);
[k0,n,x,y,GGG,parite,MMat]=deal(init{:});nby=size(MMat,1);
GGGG=zeros(length(xi),length(parite));

dx=3*mean(diff(x));fx=find(x>=min(xi)-dx&x<=max(xi)+dx);% pour gain de temps
    yyi=retelimine(yi);
    if length(yyi)==1 & abs(yyi(1))<y(nby+1);
    G=reshape(permute(GGG(fx,2:nby+1,:),[1,3,2]),[],nby)*((cos(2*(0:nby-1)*acos(yyi/y(nby+1)))*MMat)).';
    G=reshape(G,length(fx),[]);
    for ipq=1:length(parite);GGGG(:,ipq)=interp1(x(fx),G(:,ipq),xi,'spleen');if parite(ipq)==1;GGGG(:,ipq)=GGGG(:,ipq).*yyi;end;end;
    %for ipq=1:length(parite);GGGG(:,ipq)=interp1(x(fx),real(G(:,ipq)),xi,'pchip')+i*interp1(x(fx),imag(G(:,ipq)),xi,'pchip');if parite(ipq)==1;GGGG(:,ipq)=GGGG(:,ipq).*yyi;end;end;
    else;
dy=3*mean(diff(y));fy=find(y>=min(abs(yi))-dy & y<=max(abs(yi))+dy);
for ipq=1:length(parite);GGGG(:,ipq)=interpn(x(fx),y(fy),GGG(fx,fy,ipq),xi,abs(yi),'spleen');if parite(ipq)==1;GGGG(:,ipq)=GGGG(:,ipq).*yi;end;end;
     end;
     
e=zeros(length(xi),3,size(S,2));
mu=k0*(n^pol);
for i_S=1:size(S,2);
e(:,1,i_S)=-mu*S(1,i_S)*GGGG(:,1)-S(3,i_S)*GGGG(:,2)+S(2,i_S)*GGGG(:,3);
e(:,2,i_S)=-S(1,i_S)*GGGG(:,3)+(1/mu)*(-S(3,i_S)*GGGG(:,6)+S(2,i_S)*GGGG(:,5));
e(:,3,i_S)=S(1,i_S)*GGGG(:,2)+(1/mu)*(S(3,i_S)*GGGG(:,4)-S(2,i_S)*GGGG(:,6));
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function e=retpoint_periodique_conique(ep,k0,gama,S,centre,x,y,z,d,beta0,zmax,clone,parm);

Nrab=1;% diracs de chaque cote
if nargin<12;parm=[];end;if nargin<11;clone=0;end;if nargin<10;zmax=0;end;
defaultopt=struct('point',1);
point=retoptimget(parm,'point',defaultopt,'fast');
n=retsqrt(ep/k0,-1);
zmax=max([max(abs(z(:)-centre(3))),zmax,1.e-1*abs(real(2*pi/(k0*n)))]);
init=retpoint_periodique_conique_store(d,n,beta0,gama,k0,zmax,Nrab);

xmesh=imag(clone)~=0;clone=real(clone);
if ~xmesh;nx=length(x);ny=length(y);nz=length(z);[z,x,y]=ndgrid(z,x,y);end;
x=x(:)-centre(1);y=y(:)-centre(2);z=z(:)-centre(3);
p1=round(x/(2*d(1)));p2=round(y/(2*d(2)));xx=x-2*d(1)*p1;yy=y-2*d(2)*p2;% periodisation
PP=retelimine(p1+1i*p2);

if min(size(S))==1;S=S(:);end;nb_S=size(S,2);% S classe par colonnes
if clone==0;S(1:3,:)=i*S(1:3,:);end;% clonage de la source
e=retpoint_periodique_conique_interpol(init,xx,yy,z,S);

for ii=1:length(PP);f=find((p1+1i*p2)==PP(ii));if ~isempty(f);
if PP(ii)~=0;
offset_x=2*real(PP(ii))*d(1);offset_y=2*imag(PP(ii))*d(2);amp=exp(1i*beta0(1)*offset_x+1i*beta0(2)*offset_y);
for ix=-Nrab:Nrab;for iy=-Nrab:Nrab;% on rajoute les sources laterales et le centre dans les coordonnees shiftÈes (xx(f) ...
e(f,:,:)=e(f,:,:)+exp(1i*beta0(1)*ix*d(1)+1i*iy*beta0(2)*d(2))*retpoint(ep,k0,gama,S,[ix*d(1),iy*d(2),0],xx(f),yy(f),z(f),1+i);      
end;end;
e(f,:,:)=amp*e(f,:,:)-retpoint(ep,k0,gama,S,[0,0,0],x(f),y(f),z(f),1+i);% on enleve le centre dans les coordonnees non shiftÈes
else;
for ix=-Nrab:Nrab;for iy=-Nrab:Nrab;% on rajoute les sources laterales et par centre
if ix~=0 | iy~=0;e(f,:,:)=e(f,:,:)+exp(1i*beta0(1)*ix*d(1)+1i*iy*beta0(2)*d(2))*retpoint(ep,k0,gama,S,[ix*d(1),iy*d(2),0],x(f),y(f),z(f),1+i);end;      
end;end;
end;end;end;% ii

%     for ix=-Nrab:Nrab;for iy=-Nrab:Nrab;% on rajoute les sources laterales
%     %if ix~=0 | iy~=0;for ii=1:nb_S;e(:,:,ii)=e(:,:,ii)+exp(1i*beta0(1)*ix*d(1)+1i*iy*beta0(2)*d(2))*retpoint(ep,k0,gama,S(:,ii),[ix*d(1),iy*d(2),centre(3)],xx,yy,z,1+i);end;end;        
%     if ix~=0 | iy~=0;e=e+exp(1i*beta0(1)*ix*d(1)+1i*iy*beta0(2)*d(2))*retpoint(ep,k0,gama,S,[ix*d(1),iy*d(2),0],x,y,z,1+i);end;       
%     end;end;

if point==1;e=e+retpoint(ep,k0,gama,S,[0,0,0],x,y,z,1+i);end;% on ajoute la source non periodisÈe % modif 11 2012

% [offset_x,prv,ii]=retelimine(x(:)-centre(1)-xx(:));[offset_y,prv,jj]=retelimine(y(:)-centre(2)-yy(:));% modif 4 2013 pseudo periodisation
% if max(abs(offset_x))>100*eps*d(1) | max(abs(offset_y))>100*eps*d(2);amp_x=exp(1i*beta0(1)*offset_x);amp_y=exp(1i*beta0(2)*offset_y);for kk=1:nb_S;e(:,:,kk)=retdiag(amp_x(ii).*amp_y(jj))*e(:,:,kk);end;end;

if clone==0;e(:,4:6,:)=-1i*e(:,4:6,:);end;% declonage du champ
if ~xmesh;e=reshape(e,nz,nx,ny,6,nb_S);end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function init=retpoint_periodique_conique_init(d,n,beta0,gama,k0,zmax,Nrab);

k0n=k0*n;
NN_min=64;
NNx=64;NNy=128;
nnx=ceil(NNx/NN_min);nny=ceil(NNy/NN_min);
u=beta0(1)+2*pi/d(1)*(-NNx:NNx-1);x=(d(1)/(2*NNx))*(-2*NNx:nnx:2*NNx);
v=beta0(2)+2*pi/d(2)*(-NNy:NNy-1);y=(d(2)/(2*NNy))*(-2*NNy:nny:2*NNy);
[U,V]=ndgrid(u,v);
GGG=zeros(length(y),length(6),6);

[X,Y,R2,numR,GU,GV,GW]=deal(cell(2*Nrab+1,2*Nrab+1));

for ii=-Nrab:Nrab;for jj=-Nrab:Nrab;
[X{ii+Nrab+1,jj+Nrab+1},Y{ii+Nrab+1,jj+Nrab+1}]=ndgrid(x-ii*d(1),y-jj*d(2));[R2{ii+Nrab+1,jj+Nrab+1},prv,numR{ii+Nrab+1,jj+Nrab+1}]=retelimine(X{ii+Nrab+1,jj+Nrab+1}.^2+Y{ii+Nrab+1,jj+Nrab+1}.^2,10*eps);% pour gain de temps
end;end;
for iz=1:length(z);Z=z(iz);

while nnx>1;
ermax=max([abs(u(NNx/2+1))^2*exp(-(Z*abs(u(NNx/2+1)))),abs(u(NNx/2+NNx))^2*exp(-(Z*abs(u(NNx/2+NNx))))]);   
if ermax>1.e-20;break;end;
nnx=nnx/2;NNx=NNx/2;u=u(NNx+1:NNx+2*NNx);khi=khi(NNx+1:NNx+2*NNx,:);
end
while nny>1;
ermax=max([abs(v(NNy/2+1))^2*exp(-(Z*abs(v(NNy/2+1)))),abs(v(NNy/2+NNy))^2*exp(-(Z*abs(v(NNy/2+NNy))))]);   
if ermax>1.e-20;break;end;
nny=nny/2;NNy=NNy/2;v=v(NNy+1:NNy+2*NNy);khi=khi(:,NNy+1:NNy+2*NNy);
end   
% 
%     for ii=-Nrab:Nrab;for jj=-Nrab:Nrab;
%     RR=sqrt(R2{ii+Nrab+1,jj+Nrab+1}+Z^2);
%     GU{ii+Nrab+1,jj+Nrab+1}=exp(1i*k0n*RR)./(-4*pi*RR);
%     GV{ii+Nrab+1,jj+Nrab+1}=(1i*k0n./RR-1./(RR.^2)).*GU{ii+Nrab+1,jj+Nrab+1};
%     GW{ii+Nrab+1,jj+Nrab+1}=(-k0n^2./(RR.^2)-3i*k0n./(RR.^3)+3./(RR.^4)).*GU{ii+Nrab+1,jj+Nrab+1};
%     GU{ii+Nrab+1,jj+Nrab+1}=reshape(GU{ii+Nrab+1,jj+Nrab+1}(numR{ii+Nrab+1,jj+Nrab+1}),size(X{ii+Nrab+1,jj+Nrab+1}));GV{ii+Nrab+1,jj+Nrab+1}=reshape(GV{ii+Nrab+1,jj+Nrab+1}(numR{ii+Nrab+1,jj+Nrab+1}),size(X{ii+Nrab+1,jj+Nrab+1}));GW{ii+Nrab+1,jj+Nrab+1}=reshape(GW{ii+Nrab+1,jj+Nrab+1}(numR{ii+Nrab+1,jj+Nrab+1}),size(X{ii+Nrab+1,jj+Nrab+1}));
%     end;end;
% 

FF=(1/(2i*prod(d)))*exp(1i*khi*abs(Z));
FFF=FF./khi;
% disp(rettexte(NNx,NNy,nnx,nny));
for ipqr=1:size(pqr,1);
switch ipqr;
case 1;F=FFF;	
case 2;F=1i*(retdiag(u)*FFF);
case 3;F=1i*(FFF*retdiag(v));
case 4;F=1i*FF;
case 5;F=-(retdiag(u.^2)*FFF);
case 6;F=-(FFF*retdiag(v.^2));
case 7;F=-FF.*khi;
case 8;F=-retdiag(u)*FFF*retdiag(v);
case 9;F=-(retdiag(u)*FF);
case 10;F=-(FF*retdiag(v));
end;
F=(4*NNx*NNy)*(ifft2(ifftshift(F)));
GGG(:,:,iz,ipqr)=K.*F([1:nnx:2*NNx,1:nnx:2*NNx,1],[1:nny:2*NNy,1:nny:2*NNy,1]);clear F;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for ii=-Nrab:Nrab;for jj=-Nrab:Nrab;A=exp(1i*beta0(1)*ii*d(1)+1i*beta0(2)*jj*d(2));
switch ipqr;
case 1;GGG(:,:,iz,ipqr)=GGG(:,:,iz,ipqr)-A*GU{ii+Nrab+1,jj+Nrab+1};	
case 2;GGG(:,:,iz,ipqr)=GGG(:,:,iz,ipqr)-A*X{ii+Nrab+1,jj+Nrab+1}.*GV{ii+Nrab+1,jj+Nrab+1};
case 3;GGG(:,:,iz,ipqr)=GGG(:,:,iz,ipqr)-A*Y{ii+Nrab+1,jj+Nrab+1}.*GV{ii+Nrab+1,jj+Nrab+1};
case 4;GGG(:,:,iz,ipqr)=GGG(:,:,iz,ipqr)-A*Z.*GV{ii+Nrab+1,jj+Nrab+1};
case 5;GGG(:,:,iz,ipqr)=GGG(:,:,iz,ipqr)-A*GV{ii+Nrab+1,jj+Nrab+1}-A*X{ii+Nrab+1,jj+Nrab+1}.^2.*GW{ii+Nrab+1,jj+Nrab+1};
case 6;GGG(:,:,iz,ipqr)=GGG(:,:,iz,ipqr)-A*GV{ii+Nrab+1,jj+Nrab+1}-A*Y{ii+Nrab+1,jj+Nrab+1}.^2.*GW{ii+Nrab+1,jj+Nrab+1};
case 7;GGG(:,:,iz,ipqr)=GGG(:,:,iz,ipqr)-A*GV{ii+Nrab+1,jj+Nrab+1}-A*Z.^2.*GW{ii+Nrab+1,jj+Nrab+1};
case 8;GGG(:,:,iz,ipqr)=GGG(:,:,iz,ipqr)-A*X{ii+Nrab+1,jj+Nrab+1}.*Y{ii+Nrab+1,jj+Nrab+1}.*GW{ii+Nrab+1,jj+Nrab+1};
case 9;GGG(:,:,iz,ipqr)=GGG(:,:,iz,ipqr)-A*X{ii+Nrab+1,jj+Nrab+1}.*Z.*GW{ii+Nrab+1,jj+Nrab+1};
case 10;GGG(:,:,iz,ipqr)=GGG(:,:,iz,ipqr)-A*Y{ii+Nrab+1,jj+Nrab+1}.*Z.*GW{ii+Nrab+1,jj+Nrab+1};
end;
end;end;% ii jj
end;end;% ipqr iz
clear F FF FFF GU GV GW U V khi;

%%% extrapolation au centre
G0=zeros(length(x),length(y),1,10);
% a1=z(2)^2/(z(2)^2-z(1)^2);a2=-z(1)^2/(z(2)^2-z(1)^2);
% G0(:,:,:,[1,2,3,5,6,7,8])=a1*GGG(:,:,1,[1,2,3,5,6,7,8])+a2*GGG(:,:,2,[1,2,3,5,6,7,8]);% les autres composantes sont nulles

%deg=min(5,length(z));Mat=repmat(z(:).^2,1,deg).^repmat(0:deg-1,length(z),1)\eye(length(z));
%Mat=zeros(length(z));for jj=1:length(z);Mat(:,jj)=cos(2*(jj-1)*acos(z.'/z(end)));end;Mat=(-1).^[0:length(z)-1]/Mat;% polynomes de Tchebitchev
%MMat=zeros(length(z));for jj=1:length(z);MMat(:,jj)=cos(2*(jj-1)*acos(z.'/z(end)));end;MMat=inv(MMat);Mat=(-1).^[0:length(z)-1]*MMat;% polynomes de Tchebitchev
%prvMat=Mat;prvMMat=MMat;

MMat=zeros(nbz);for jj=1:nbz;MMat(:,jj)=cos(2*(jj-1)*acos(z(1:nbz).'/z(nbz)));end;MMat=inv(MMat);Mat=(-1).^[0:nbz-1]*MMat;% polynomes de Tchebitchev

%for ipqr=[1,2,3,5,6,7,8];for ii=1:length(z);G0(:,:,1,ipqr)=G0(:,:,1,ipqr)+GGG(:,:,ii,ipqr)*Mat(1,ii);end;end;

GGG=reshape(GGG,length(x)*length(y),length(z),10);
for ipqr=[4,9,10];GGG(:,:,ipqr)=GGG(:,:,ipqr)*retdiag(1./z);end;% on divise par z les impairs pour les rendre pairs
GGG=reshape(GGG,length(x),length(y),length(z),10);
%prvG0=G0;for ipqr=1:10;for ii=1:length(z);prvG0(:,:,1,ipqr)=prvG0(:,:,1,ipqr)+GGG(:,:,ii,ipqr)*prvMat(1,ii);end;end;
for ipqr=1:10;for ii=1:nbz;G0(:,:,1,ipqr)=G0(:,:,1,ipqr)+GGG(:,:,ii,ipqr)*Mat(1,ii);end;end;

% p=20;t=linspace(-1,1,2*p);t=t(p+1:2*p);Tch=zeros(p);for jj=1:p;Tch(:,jj)=cos(2*(jj-1)*acos(t.'));end;M=Tch\eye(p);
% symetrisation en z
parite=mod(pqr(:,end),2);

%for ii=1:10;fx=find(abs(x)<max(x)/4);fy=find(abs(y)<max(y)/4);retfig;figure;retcolor(y(fy),x(fx),real(GGG(fx,fy,1,ii)),2i);hold on;;pause;end;stop


N=2;Ix0=retminabs(x);Iy0=retminabs(y);fx=Ix0-N:Ix0+N;fy=Iy0-N:Iy0+N;
[G00,G10,G01,G20,G11,G02]=deal(zeros(nbz,length(parite)));
for ipqr=1:length(parite);for iz=1:nbz; 
GGGx=retderivee(x(fx),GGG(fx,fy,iz,ipqr),1);GGGy=retderivee(y(fy),GGG(fx,fy,iz,ipqr),2);
GGGxx=retderivee(x(fx),GGGx,1);GGGxy=.5*(retderivee(x(fx),GGGy,1)+retderivee(y(fy),GGGx,2));GGGyy=retderivee(y(fy),GGGy,2);
G00(iz,ipqr)=GGG(Ix0,Iy0,iz,ipqr);G10(iz,ipqr)=GGGx(N+1,N+1);G01(iz,ipqr)=GGGy(N+1,N+1);G11(iz,ipqr)=GGGxy(N+1,N+1);G20(iz,ipqr)=GGGxx(N+1,N+1);G02(iz,ipqr)=GGGyy(N+1,N+1);
end;end;
z=[0,z];
GGG=cat(3,G0,GGG);

init={k0,n,x,y,z,GGG,parite,MMat,G00,G10,G01,G20,G11,G02};

