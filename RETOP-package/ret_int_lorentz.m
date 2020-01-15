function f=ret_int_lorentz(varargin);
% calcul de: somme(E1^H2-E2^H1)
%         1 D            
%**************************
% f=ret_int_lorentz(e1,e2,wx,wy);
% e1,e2,wx,wy  calcules par  [e,y,wy]=retchamp(init,a,sh,sb,inc,x,...				
% wx poids pour integration de Gauss en x
% wy poids pour integration de Gauss en y
%
%        2 D           
% **************************
% f=ret_int_lorentz(e1,e2,wx,wy,wz); 
% 
% e1,e2,wx,wy,wz  calcules par  [e,z,wz]=retchamp(init,a,sh,sb,inc,x,...				
% wx poids pour integration de Gauss en x
% wy poids pour integration de Gauss en y
% wz poids pour integration de Gauss en z
%
% Dans les 2 cas, UN DES POIDS EST VIDE POUR INDIQUER QUE L'ON INTEGRE PERPENDICULAIREMENT A CETTE DIRECTION
%                  si les poids sont absents en entrée f=E1^H2-E2^H1



% See also: RETLORENTZ RETPOYNTING
sz=size(varargin{1});
if sz(end)<6;f=ret_int_lorentz_1D(varargin{:});else;f=ret_int_lorentz_2D(varargin{:});end;
%if length(size(varargin{1}))<4;f=ret_int_lorentz_1D(varargin{:});else;f=ret_int_lorentz_2D(varargin{:});end;
%if nargin<5;f=ret_int_lorentz_1D(varargin{:});else;f=ret_int_lorentz_2D(varargin{:});end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function f=ret_int_lorentz_1D(e1,e2,wx,wy);
if nargin<3;
if size(e1,3)>2;
f=cat(3,-e1(:,:,1).*e2(:,:,3)+e2(:,:,1).*e1(:,:,3) ,    e1(:,:,1).*e2(:,:,2)-e2(:,:,1).*e1(:,:,2));
else;
f=[-e1(:,1).*e2(:,3)+e2(:,1).*e1(:,3) ,    e1(:,1).*e2(:,2)-e2(:,1).*e1(:,2)];
end;
return;
end;

if isempty(wy);% integration perpendiculaire à  oy % - E1z * H2y + E2z * H1y
f=((e1(:,:,1).*e2(:,:,2)-e2(:,:,1).*e1(:,:,2))*wx(:)).';%  E1z * H2x - E2z * H1x
else;          % integration perpendiculaire à  ox 
f=wy(:).'*(-e1(:,:,1).*e2(:,:,3)+e2(:,:,1).*e1(:,:,3));% - E1z * H2y + E2z * H1y
end;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function f=ret_int_lorentz_2D(e1,e2,wx,wy,wz);
sz=size(e1);
if nargin<3;
if size(e1,4)>4;
    f=cat(4,e1(:,:,:,2).*e2(:,:,:,6)-e1(:,:,:,3).*e2(:,:,:,5)-e2(:,:,:,2).*e1(:,:,:,6)+e2(:,:,:,3).*e1(:,:,:,5),...
    e1(:,:,:,3).*e2(:,:,:,4)-e1(:,:,:,1).*e2(:,:,:,6)-e2(:,:,:,3).*e1(:,:,:,4)+e2(:,:,:,1).*e1(:,:,:,6),...
    e1(:,:,:,1).*e2(:,:,:,5)-e1(:,:,:,2).*e2(:,:,:,4)-e2(:,:,:,1).*e1(:,:,:,5)+e2(:,:,:,2).*e1(:,:,:,4));
else;
%     f=[e1(:,:,:,2).*e2(:,:,:,6)-e1(:,:,:,3).*e2(:,:,:,5)-e2(:,:,:,2).*e1(:,:,:,6)+e2(:,:,:,3).*e1(:,:,:,5),...
%     e1(:,:,:,3).*e2(:,:,:,4)-e1(:,:,:,1).*e2(:,:,:,6)-e2(:,:,:,3).*e1(:,:,:,4)+e2(:,:,:,1).*e1(:,:,:,6),...
%     e1(:,:,:,1).*e2(:,:,:,5)-e1(:,:,:,2).*e2(:,:,:,4)-e2(:,:,:,1).*e1(:,:,:,5)+e2(:,:,:,2).*e1(:,:,:,4)];

%     f=[e1(:,2).*e2(:,6)-e1(:,3).*e2(:,5)-e2(:,2).*e1(:,6)+e2(:,3).*e1(:,5),...
%     e1(:,3).*e2(:,4)-e1(:,1).*e2(:,6)-e2(:,3).*e1(:,4)+e2(:,1).*e1(:,6),...
%     e1(:,1).*e2(:,5)-e1(:,2).*e2(:,4)-e2(:,1).*e1(:,5)+e2(:,2).*e1(:,4)];
     f=zeros(size(e1,1),3);% gain de temps
     f(:,1)=e1(:,2).*e2(:,6)-e1(:,3).*e2(:,5)-e2(:,2).*e1(:,6)+e2(:,3).*e1(:,5);
     f(:,2)=e1(:,3).*e2(:,4)-e1(:,1).*e2(:,6)-e2(:,3).*e1(:,4)+e2(:,1).*e1(:,6);
     f(:,3)=e1(:,1).*e2(:,5)-e1(:,2).*e2(:,4)-e2(:,1).*e1(:,5)+e2(:,2).*e1(:,4);;
end;
return;
end;

if isempty(wx);% integration perpendiculaire à ox 
% f=zeros(sz(2),1);prv=reshape(wz(:)*wy(:).',sz(1),1,sz(3));	
% for ii=1:sz(2);% +E1y * H2z -E1z * H2y -E2y * H1z + E2z * H1y
% f(ii)=sum(retcolonne(prv.*(e1(:,ii,:,2).*e2(:,ii,:,6)-e1(:,ii,:,3).*e2(:,ii,:,5)-e2(:,ii,:,2).*e1(:,ii,:,6)+e2(:,ii,:,3).*e1(:,ii,:,5))));
% end;
% gain de temps modif 3 2014
f=reshape(permute((e1(:,:,:,2).*e2(:,:,:,6)-e1(:,:,:,3).*e2(:,:,:,5)-e2(:,:,:,2).*e1(:,:,:,6)+e2(:,:,:,3).*e1(:,:,:,5)),[2,1,3]),sz(2),sz(1)*sz(3))*retcolonne(wz(:)*wy(:).');
end;

if isempty(wy);% integration perpendiculaire à oy 
% f=zeros(sz(3),1);prv=reshape(wz(:)*wx(:).',sz(1),sz(2),1);	
% for ii=1:sz(3);% +E1z * H2x -E1x * H2z -E2z * H1x + E2x * H1z
% f(ii)=sum(retcolonne(prv.*(e1(:,:,ii,3).*e2(:,:,ii,4)-e1(:,:,ii,1).*e2(:,:,ii,6)-e2(:,:,ii,3).*e1(:,:,ii,4)+e2(:,:,ii,1).*e1(:,:,ii,6))));
% end;
f=reshape(permute((e1(:,:,:,3).*e2(:,:,:,4)-e1(:,:,:,1).*e2(:,:,:,6)-e2(:,:,:,3).*e1(:,:,:,4)+e2(:,:,:,1).*e1(:,:,:,6)),[3,1,2]),sz(3),sz(1)*sz(2))*retcolonne(wz(:)*wx(:).');

end;
if isempty(wz);% integration perpendiculaire à oz
% f=zeros(sz(1),1);prv=reshape(wx(:)*wy(:).',1,sz(2),sz(3));
% for ii=1:sz(1);% +E1x * H2y -E1y * H2x -E2x * H1y + E2y * H1x
% f(ii)=sum(retcolonne(prv.*(e1(ii,:,:,1).*e2(ii,:,:,5)-e1(ii,:,:,2).*e2(ii,:,:,4)-e2(ii,:,:,1).*e1(ii,:,:,5)+e2(ii,:,:,2).*e1(ii,:,:,4))));
% end;
f=reshape((e1(:,:,:,1).*e2(:,:,:,5)-e1(:,:,:,2).*e2(:,:,:,4)-e2(:,:,:,1).*e1(:,:,:,5)+e2(:,:,:,2).*e1(:,:,:,4)),sz(1),sz(2)*sz(3))*retcolonne(wx(:)*wy(:).');


end;
