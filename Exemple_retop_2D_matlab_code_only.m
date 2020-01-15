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


%%%%%%% Exemple_retop_2D_matlab_code_only.m
%%%%%%% This is a Matlab script only, which executes three major tasks :
%%%%%%% 1) Solve the Maxwell's equation for a dipole line souce (2D problem) emitting in a uniform medium (use the retpoint.m function)  
%%%%%%% 2) calculate the points on a box around the dipole and save their coordinates in prv_coordonnees.mat 
%%%%%%% 3) calculate the electromagnetic fields at the points on the box and send the fields into RETOP for retrieving the far-field Radiation Diagram

% specific parameters (user does not change them)  
option_cal_champ=2;option_i=1;

% physical parameters
ld=.6e-6; % wavelength in m
k0=2*pi/ld;
n_strates=1.5;z_strates=[];nh=n_strates(1); % define the stratified medium (uniform medium of refractive index 1.5 in the example). See documentation.
pol=2; % polarisation 0: E//          2: H//

% definition of the box in COMSOL coordinates (the second coordinate, y, is perpendicular to the layers)
Centre=[0,0];L=[3*ld,2*ld];N_gauss=[150;140];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Thanks to the function extract_pointcoordinate_on_box.m,
% we compute the coordinates on the box, order them and save them in prv_coordonnees.mat
init=retop(@extract_pointcoordinate_on_box,pol,k0,n_strates,z_strates,Centre,L,N_gauss,struct('option_cal_champ',option_cal_champ,'option_i',option_i));


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% compute the field at the points defined by the coordinates [x,y]=[coordonnees{1} coordonnees{2}]
% the exemple compute the field radiated by the dipole line source S_inc located at C_inc, and use a matlab function retpoint_comsol  
load prv_coordonnees


% S_inc  dipole line source (3 component vector)
% for pol=0, S_inc(1)=pz,S_inc(2)=µx,S_inc(3)=µy
% for pol=2, S_inc(1)=µz,S_inc(2)=px,S_inc(3)=py)
S_inc=[0 1 0]; % µx or px line source for TE or TM polarization
C_inc=[0,0];   % x,y coordinates of the dipole line source

% cal_champ=@(x,y) ret2comsol(retpoint(k0*1.e-6*n_strates(1)^2,k0*1.e-6,S_inc,C_inc,x*1.e6,y*1.e6,pol,1i),option_i,pol); % compute the field radiated by the dipole S_inc
% champ=cal_champ(coordonnees{:});                                                          % compute the field at the coordinates on the box

% cal_champ=@(x,y,z) ret2comsol(retpoint(k0*1.e-6*n_strates(1)^2,k0*1.e-6,S_inc,C_inc,x*1.e6,y*1.e6,z*1.e6,1i),option_i,1); % transform the field calculated with retpoint in SI
% champ=cal_champ(coordonnees{:});                                              % compute the field at the coordinates on the box

champ=retpoint_comsol(k0*n_strates(1)^2,k0,S_inc,C_inc,coordonnees{1},coordonnees{2},pol);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% compute the plane wave decomposition in free space and make some tests
init=retop(@(x,y) champ,pol,k0,n_strates,z_strates,Centre,L,N_gauss,struct('option_cal_champ',option_cal_champ,'option_i',option_i));

teta=retgauss(-pi/2,pi/2,200);
u=sin(teta);v=cos(teta);
uh=k0*n_strates(1)*u;ub=uh;
angles_h=retop(init,uh,1,struct('test',1));% plane-wave decomposition above
angles_b=retop(init,ub,-1);% plane-wave decomposition below

figure;hold on;
plot(u.*angles_h.F.',v.*angles_h.F.',u.*angles_b.F.',-v.*angles_b.F.','-k','linewidth',2);
axis equal;xlabel('x');ylabel('y');title('radiation diagram in the x-y plane')

