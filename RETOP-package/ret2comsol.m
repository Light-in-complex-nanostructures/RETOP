function e=ret2comsol(e,option_SI,pol);% e ret--> e comsol
if option_SI>0;e=conj(e);end;
ratioH=19.409541814833513;%ratioH=sqrt(retconstantes('Z0'));
ratioE=51.5210513231055e-3;%ratioE=1./sqrt(retconstantes('Z0'));
switch(pol);
case 1;e(:,1:3)=e(:,1:3)/ratioE;e(:,4:6)=e(:,4:6)/ratioH;% 3D 
case 0;e(:,1)=e(:,1)/ratioE;e(:,2:3)=e(:,2:3)/ratioH;% TE 
case 2;e(:,1)=e(:,1)/ratioH;e(:,2:3)=-e(:,2:3)/ratioE;% TM 
end;
