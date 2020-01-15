function font=retfont(k,type);
% font=retfont(k,type);
%  font=retfont ==> font={'interpreter','none','fontname','Times NewRoman','fontweigh','bold'}
% d'après le programme a_plot de Antony on peut tout appliquer :retfont(gca), ou retfont(gca,font)
% trace avec des types de lignes engendrés automatiquement
%  si type=1 lignes  '-k',':r','-.g'
%  si type=2 lignes avec marqueurs 'ok','xr',...
%  si type=3 lignes avec marqueurs '.-k','o:r','x-.g',
% 
% appliquer 'font={'interpreter','none','fontname','Times New Roman','fontweigh','bold'}; à tous les subplots d'une figure: retfont(gcf,0);
% font=[retfont,{'fontsize',15}],appliquer 'font' à tous les subplots d'une figure: retfont(gcf,font);
% appliquer 'font standart ' à un axis: retfont(gca);
% appliquer 'font' à un axis: retfont(gca,font);
%
%%% EXEMPLE DE TRACE
% figure;
% subplot(3,1,1);hold on;for ii=1:5;plot(rand(1,5),retfont(ii));end;retfont(gca);
% subplot(3,1,2);hold on;for ii=1:5;plot(rand(1,5),retfont(ii,2));end;retfont(gca);
% subplot(3,1,3);leg={};hold on;for ii=1:5;plot(rand(1,5),retfont(ii,3));leg=[leg,{['ii=',num2str(ii)]}];end;legend(leg{:});
% retfont(gcf,0);




if nargin<1;font={'interpreter','none','fontname','Times New Roman','fontweigh','bold'};return;end;

if nargin<2;type=1;end;


if isnumeric(k);kk=k;else;switch get(k(1),'type');case 'figure';kk=0;otherwise kk=.5;end;end;% compatibilite 2014
if  abs(kk(1)-round(kk(1)))>100*eps;% k est le handel d'un plot
if nargin<2;type=retfont;end;% d'après le programme a_plot de Antony 2011
set(get(k,'Title'),type{:});
set(get(k,'Xlabel'),type{:});
set(get(k,'Ylabel'),type{:});
set(get(k,'Zlabel'),type{:});
set(k,type{3:end});
return;
else
if iscell(type);ch=get(k,'children');for ii=1:length(ch);switch get(ch(ii),'type');case 'axes';retfont(ch(ii),type);end;end;return;end
if type==0;ch=get(k,'children');for ii=1:length(ch);switch get(ch(ii),'type');case 'axes';retfont(ch(ii));end;end;return;end
% return
% 	
% end;

% lin={{'-k','-r','-g','-b','-m','-c','--k','--r','--g','--b','--m','--c','-.k','-.r','-.g','-.b','-.m','-.c',':k',':r',':g',':b',':m',':c'};...
% 	{'.k','.r','.g','.b','.m','.c','*k','*r','*g','*b','*m','*c','+k','+r','+g','+b','+m','+c','xk','xr','xg','xb','xm','xc'}};
lin={{'-k',':r','-.g','--b','-m',':c','-.k','--r','-g',':b','-.m','--c','-k',':r','-.g','--b','-m',':c','-.k','--r','-g',':b','-.m','--c'};...
	 {'ok','xr','+g','*b','sm','dc','vk','^r','<g','>b','pm','hc','ok','xr','+g','*b','sm','dc','vk','^r','<g','>b'};...
	{'.-k','o:r','x-.g','+--b','*-m','s:c','d-.k','v--r','^-g','<:b','>-.m','p--c','h-k','.:r','o-.g','x--b','*-m','s:c','d-.k','v--r','^-g','<:b','>-.m','p--c'}};
	
font=lin{type}{mod(k,length(lin{type})-1)+1};
return
	
end;


