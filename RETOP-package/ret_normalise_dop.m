function ret_normalise_dop(Fac);

% Il arrive qu'une figure 3D tracée avec 'surface' apparaisse en noir. Ceci est du à des valeurs trop grandes (1.e20) ou trop petites (1.e-20)
% ret_normalise_dop(Fac) aprés selection d'un axis, permet le multiplier par Fac 
% Attention, les Xtick restent inchangés mis à part une puissance de 10 que je ne sais pas controler ..
% 
% ret_normalise_dop(Fac);
% renormalisation d'un DOP trop grand ou trop petit Multiplication par Fac
%set(gca,'XTicklabelmode','manual');set(gca,'YTicklabelmode','manual');set(gca,'ZTicklabelmode','manual');

set(gca,'Xlim',get(gca,'Xlim')*Fac,'XTick',get(gca,'XTick')*Fac,'XTickLabel',get(gca,'XTickLabel'),...
'Ylim',get(gca,'Ylim')*Fac,'YTick',get(gca,'YTick')*Fac,'YTickLabel',get(gca,'YTickLabel'),...
'Zlim',get(gca,'Zlim')*Fac,'ZTick',get(gca,'ZTick')*Fac,'ZTickLabel',get(gca,'ZTickLabel'));

ch=get(gca,'children');
for ii=1:length(ch);
switch get(ch(ii),'type');case 'surface'
set(ch(ii),'Cdata',get(ch(ii),'Cdata')*Fac);
set(ch(ii),'Xdata',get(ch(ii),'Xdata')*Fac);
set(ch(ii),'Ydata',get(ch(ii),'Ydata')*Fac);
set(ch(ii),'Zdata',get(ch(ii),'Zdata')*Fac); 
end;end;

%set(gca,'XTicklabelmode','manual');set(gca,'YTicklabelmode','manual');set(gca,'ZTicklabelmode','manual');



