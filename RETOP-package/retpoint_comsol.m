function e=retpoint_comsol(ep,mu,s,xs,x,y,z,clone);
if nargin<8;clone=1i;end;
ep=1.e-6*ep;mu=1.e-6*mu;x=1.e6*x;y=1.e6*y;% unites reticolo (m--> micron)
if length(xs)==3;z=1.e6*z;pol=1;else;pol=z;end% 3D en 2D z=pol
option_i=1;
e=retpoint(ep,mu,s,xs,x,y,z,clone);
e=ret2comsol(e,option_i,pol);
