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

%%%%%%% This is a main Matlab script, which runs in the Matlab-COMSOL livelink.
%%%%%%% It executes three major tasks :
%%%%%%% 1) Solve the Maxwell's equation of the light-emission problem
%%%%%%% 2) Send the electromagnetic field on a 'box' into RETOP for retrieving the Free-space Radiation Diagram
%%%%%%% 3) Send the electromagnetic field on a 'box' into RETOP for retrieving the Guided-mode Radiation Diagram

clear, close all

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Part 1: Solver the COMSOL model (in the Matlab-COMSOL enviroment)
disp('the computation lasts about 100 seconds')

% initialize the COMSOL model
import com.comsol.model.*
import com.comsol.model.util.*
ModelUtil.showProgress(true)% show the panel of progression
% load the model, where the 'box' is defined with a size (Lx,Ly,Lz)
nom='double_dipole_in_waveguide_web.mph';  % the name of the COMSOL model sheet
model = mphload(nom);  
% Solving the COMSOL model
% wavelength=1e-6; model.param.set('lambda',num2str(wavelength)); % change the wavelength in COMSOL calculation (if needed)
model.geom('geom1').run(); % generate the geometry
model.mesh('mesh1').run(); % generate the mesh
model.study('std1').run(); % solve the Maxwell's equation in the model 
%-------> the COMSOL model has been solved here


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Part 2: Retrieve the Free-Space radiation pattern
% the input and output have to be in SI unit (like in the commercial software COMSOL);
tic;

%%% <<cal_champ>> is based on the subroutine <<extract_comsol_field>> that interfaces COMSOL with Matlab
%%% <<extract_comsol_field>> takes the EM field from the solved COMSOL model (done in Part.1)
opt_field=1;                                % Take the "total field" or "scattered field" from COMSOL. "1": total field; "0": scattered field (avaialable only in "scattering formulation"). 
cal_champ=@(x,y,z) extract_comsol_field(model,x,y,z,opt_field);
wavelength=1000e-9;                                     % "wavelength" (unit:meter)
wavenumber=2*pi/wavelength;
n_air=1;n_glass=1.5;n_sub=1.2;                  % refractive indices of the three layers (upper cladding, core, and low cladding)
refractive_indices=[n_air,n_glass,n_sub]; % the refractive indices of the stratified medium
                                                                            % n_air (upper cladding)  n_glass (core)   n_sub (substrate)
z_layers=[100e-9,-100e-9];                            % z coordinate of each discontinuity (unit: meter)
center_x=0; center_y=0; center_z=0;        
box_center=[center_x,center_y,center_z];    % center of the box (in RETOP unit : meter)    
Lx=300e-9;Ly=100e-9;Lz=300e-9;
box_size=[Lx,Ly,Lz];             % Box Size (in meter)--> the box has the same size as the one defined in COMSOL model
N_gauss=[8,8;8,8;8,8];           % A parameter for generating the sampling locations on the "box" (normally, the present setting is enough)
option_i=1;                               % 'option_i=1' : input field has convention exp(i*omega*t)
                                          % 'option_i=-1': input field has convention exp(-i*omega*t)
%%% Initialize RETOP by inputting material and geometrical parameters
init=retop(cal_champ,wavenumber,refractive_indices,z_layers,box_center,box_size,N_gauss,struct('option_cal_champ',2,'option_i',option_i));
% define the (theta --> polar angle, phi --> azimuthal angle) (for details, see the user guide)
[teta,wteta]=retgauss(0,pi/2,10,3);         % polar angle
[phi,wphi]=retgauss(0,2*pi,10,7);           % azimuthal angle
[Teta,Phi]=ndgrid(teta,phi);                % generate the grid
u=sin(Teta).*cos(Phi);                      % weight in x-coordinate
v=sin(Teta).*sin(Phi);                      % weight in y-coordinate
w=cos(Teta);                                % weight in z-coordinate

%%% Planewave decomposition in the upper free-space (if the corresponding refractive index is real)
if imag(refractive_indices(1))==0;
    uh=wavenumber*refractive_indices(1)*u;vh=wavenumber*refractive_indices(1)*v;  %  wavenumber kx=uh, ky=vh, refractive_indices(1): refractive index of the upper free-space
    direction=1;                                          %  Attention: direction=1 (upper space, +z)
    angles_up=retop(init,uh,vh,direction,struct('test',1));%  planewave decomposition for the upper space
    % the option "struct('test',1)" allows displaying the field on the box in order to check the quality of sampling
end;

%%% Planewave decomposition in the lower free-space (if the corresponding refractive index is real)
if imag(refractive_indices(end))==0;
    ub=wavenumber*refractive_indices(end)*u;vb=wavenumber*refractive_indices(end)*v;  %  wavenumber kx=ub, ky=vb, refractive_indices(end): refractive index of the lower free-space
    direction= -1;                                        %  Attention: direction= -1 (lower space, -z)
    angles_dn=retop(init,ub,vb,direction);                 %  planewave decomposition for the lower space
end;

%%% calculate the energies radiated in upper half-space & in the substrate 
[wt,wd]=ndgrid(wteta,wphi); aef_w=wt(:).*wd(:).*sin(Teta(:)); % weight of each solid angle
energie_up=sum(angles_up.F(:).*aef_w(:));                     % the total amount of power radiated in the upper half-space
energie_down=sum(angles_dn.F(:).*aef_w(:));                   % the total amount of power radiated in the substrate

%%% plot the 3D Radiation Diagram or Radiation Pattern in Free Space
figure;hold on;
Fh=reshape(angles_up.F,size(uh));              % the flux at each solid angle, upper space ("h"-->upper)
Fb=reshape(angles_dn.F,size(ub));              % the flux at each solid angle, lower space ("b"-->lower)
if imag(refractive_indices(1))==0;  surf(u.*Fh,v.*Fh,w.*Fh,'facecolor',[.5,.5,.5],'LineStyle','none');end;
if imag(refractive_indices(end))==0;surf(u.*Fb,v.*Fb,-w.*Fb,'facecolor',[.5,.5,.5],'LineStyle','none');end;
set(gca,'projection','perspective');
axis tight;axis equal;xlabel('x');ylabel('y');
lighting PHONG;light('position',[-1,0,1],'color','r');light('position',[-1,-5,1],'color','w');light('position',[-1,-5,1],'color','y');
toc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Part 3: Retrieve the guided-mode radiation patterns
tic
neff_TE=1.24;  pol_TE=0; % specify an approximate value of the effective index ("neff_TE") and the polarization ("pol_TE") of the TE0 mode
neff_TM=1.20; pol_TM=2; % specify an approximate value of the effective index ("neff_TM") and the polarization ("pol_TM") of the TM0 mode 
LLLL=-5:5;              % set the range of azimuthal indices n=[-5,5] of the cylindrical mode, [-5,5]--> 11 sub-modes are applied, see the Paper ??? for details            
[teta_mode,wteta_mode]=retgauss(0,2*pi,400); % azimuthal angle
angles_TM=retop(init,neff_TM,pol_TM,LLLL,teta_mode); % mode decomposition of the TM0 mode
angles_TE=retop(init,neff_TE,pol_TE,LLLL,teta_mode); % mode decomposition of the TE0 mode
energie_modes_TM=sum(angles_TM.F(:).*wteta_mode(:)); % the total amount of power coupled to the TM0 mode 
energie_modes_TE=sum(angles_TE.F(:).*wteta_mode(:)); % the total amount of power coupled to the TE0 mode ;

% plot the 2D guided-mode Radiation Patterns
figure;
subplot(1,2,1);hold on; plot(angles_TM.F.*cos(teta_mode),angles_TM.F.*sin(teta_mode),'-k','linewidth',2);axis equal;
subplot(1,2,2);hold on; plot(angles_TE.F.*cos(teta_mode),angles_TE.F.*sin(teta_mode),'--r','linewidth',2);axis equal;

toc

stop
% save the data
save COMSOL_Dop_DoubleDipole_SlabGuide_lam1000   wavelength refractive_indices z_strates opt_field option_cal_champ n_air L Centre N_gauss teta wteta ...
    phi wphi Teta Phi u v w uh vh ub vb angles_up angles_dn wt wd aef_w energie_pws_h energie_box ...
    LLLL teta_mode wteta_mode neff_TM angles_TM energie_modes_TM neff_TE angles_TE energie_modes_TE

