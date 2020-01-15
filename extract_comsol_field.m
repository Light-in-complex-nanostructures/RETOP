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

function EH=extract_comsol_field(model,x,y,z_or_pol,opt_field)
%%%%% extract_comsol_field.m %%%%%
% this function extracts the (scattered or total) field from COMSOL multiphysics
% 'model' : the COMSOL model that provides the field on the box
% x,y are vectors that define the x and y coordinates of the points for which we extract the COMSOL field
% 'z_or_pol' is the coordinate vector 'z' for 3D cases, or for 2D cases (since we do not have z) 'z_or_pol' is the 'polarization': a single number '0' for TE or '2' for TM
% 2D only: TE Polarization: (Ez, Hx, Hy) or TM Polarization: (Hz, Ex, Ey)
% 'opt_field': extract the total field (if opt_field==1) or scattered field (if opt_field==0)

%% 3D case %%
if length(x)==length(z_or_pol) 
    z=z_or_pol;
    coord=[x(:),y(:),z(:)];
    EH=zeros(size(coord,1),6);
    switch opt_field
        % case=1 ; total field
        case 1; [Ex,Ey,Ez,Hx,Hy,Hz]=mphinterp(model,{'emw.Ex', 'emw.Ey','emw.Ez','emw.Hx','emw.Hy','emw.Hz'},'coord',coord.','Complexout','on');
        % case=0 ; scattered field
        case 0; [Ex,Ey,Ez,Hx,Hy,Hz]=mphinterp(model,{'emw.relEx', 'emw.relEy','emw.relEz','emw.relHx','emw.relHy','emw.relHz'},'coord',coord.','Complexout','on');
    end;
    EH=([Ex;Ey;Ez;Hx;Hy;Hz].');

%% 2D case %% 
else 
    pol=z_or_pol;
    coord=[x(:),y(:)];
    EH=zeros(size(coord,1),3);
    if pol==0 % TE Polarization:(Ez, Hx, Hy)
        switch opt_field
            % case=1 ; total field
            case 1; [Ez,Hx,Hy]=mphinterp(model,{'emw.Ez','emw.Hx','emw.Hy'},'coord',coord.','Complexout','on');
            % case=0 ; scattered field
            case 0; [Ez,Hx,Hy]=mphinterp(model,{'emw.relEz','emw.relHx','emw.relHy'},'coord',coord.','Complexout','on');
        end;
        EH=[Ez;Hx;Hy].';
        
    elseif pol==2 % TM Polarization:(Hz, Ex, Ey)
        switch opt_field
            % case=1 ; total field
            case 1; [Hz,Ex,Ey]=mphinterp(model,{'emw.Hz','emw.Ex', 'emw.Ey'},'coord',coord.','Complexout','on');
            % case=0 ; scattered field
            case 0; [Hz,Ex,Ey]=mphinterp(model,{'emw.relHz','emw.relEx', 'emw.relEy'},'coord',coord.','Complexout','on');
        end;
        EH=([Hz;Ex;Ey].'); 
    end
end

end