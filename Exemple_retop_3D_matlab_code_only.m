%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initial software: NtoFField
% Version 8.0: NtoFField
% 
% Co-authors: Jianji Yang, Jean Paul Hugonin and Philippe Lalanne
%
% Owners: Centre National de la Recherche Scientifique (CNRS) and 
%         Institut d'Optique Graduate School (IOGS)
%
% Copyright (C) 2015-2016, spread under the terms and conditions of the  
% license GPL Version 3.0
%
% See gpl-3.0.txt  or  http://www.gnu.org/licenses/gpl-3.0-standalone.html 
% for a full version of the licence 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%% Exemple_retop_3D_matlab_code_only.m
%%%%%%% This is a Matlab script only, which executes three major tasks :
%%%%%%% 1) Solve the Maxwell's equation for a point-dipole source (3D problem)emitting in a uniform medium (with retpoint.m)  
%%%%%%% 2) calculate the points on a box around the dipole and save their coordinates in prv_coordonnees.mat 
%%%%%%% 3) calculate the electromagnetic fields at the points on the box and send the fields into RETOP for retrieving the far-field Radiation Diagram

% specific parameters (user does not change them)  
option_cal_champ=2;option_i=1;


% physical parameters
ld=.6e-6; % wavelength in m
k0=2*pi/ld;
n_strates=1.5;z_strates=[];nh=n_strates(1);% define the stratified medium (uniform medium of refractive index 1.5 in the example). See documentation.

% definition of the box in COMSOL coordinates (the third coordinate, z, is perpendicular to the layers)
Centre=[0,0,0];L=[2*ld,2*ld,1*ld];N_gauss=[10,5;11,5;10,10];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Thanks to the function extract_pointcoordinate_on_box.m,
% we compute the coordinates on the box, order them and save them in prv_coordonnees.mat
init=retop(@extract_pointcoordinate_on_box,k0,n_strates,z_strates,Centre,L,N_gauss,struct('option_cal_champ',option_cal_champ,'option_i',option_i));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% compute the field at the points defined by the coordinates [x,y,z]=[coordonnees{1} coordonnees{2} coordonnees{3}]
% the exemple compute the field radiated by the dipole S_inc located at C_inc, and use a matlab function retpoint_comsol  
load prv_coordonnees;

S_inc=[0 0 1 0 0 0 0];            % dipole source pz=1 (6 component vector: S_inc(1)=px,S_inc(2)=py,S_inc(3)=pz,S_inc(4)=µx,S_inc(5)=µy,S_inc(6)=µz)
C_inc=[0,0,0];                    % x,y,z coordinates of the dipole

champ=retpoint_comsol(k0*n_strates(1)^2,k0,S_inc,C_inc,coordonnees{1},coordonnees{2},coordonnees{3}); % compute the field at the coordinates on the box 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% compute the plane wave decomposition in free space and make some tests
init=retop(@(x,y,z) champ,k0,n_strates,z_strates,Centre,L,N_gauss,struct('option_cal_champ',option_cal_champ,'option_i',option_i));% calcul du champ sur la boite
teta=linspace(0,pi/2,20);phi=linspace(0,2*pi,60);
% if one fixes phi, one may get the polar plot in the different planes XY, XZ, ZY (2D diagram pattern)
[Teta,Phi]=ndgrid(teta,phi);
u=sin(Teta).*cos(Phi);v=sin(Teta).*sin(Phi);w=cos(Teta);

uh=k0*n_strates(1)*u;vh=k0*n_strates(1)*v;
angles_h=retop(init,uh,vh,1,struct('test',1));% calcul du DOP en haut
ub=uh;vb=vh;angles_b=retop(init,ub,vb,-1);    % calcul du DOP en bas

figure;hold on;% tracé du DOP
Fh=reshape(angles_h.F,size(uh));surf(u.*Fh,v.*Fh,w.*Fh,'facecolor',[.5,.5,.5],'LineStyle','none');
Fb=reshape(angles_b.F,size(ub));surf(u.*Fb,v.*Fb,-w.*Fb,'facecolor',[.5,.5,.5],'LineStyle','none');
axis tight;axis equal;xlabel('x');ylabel('y');zlabel('z');view([-20,40]);
lighting gouraud;light('position',[-1,0,1],'color','r');light('position',[-1,-5,1],'color','w');light('position',[-1,-5,1],'color','y');
title('radiation diagram');
ret_normalise_dop(1.e13)
  
