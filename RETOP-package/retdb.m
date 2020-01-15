function db=retdb(n,ld);
% db=retdb(n,ld);
% attenuation en db/mm d'un indice complexe n à la longueur d'onde ld (microns)
% 
% conversion ld MICRONS , en Tera Hertz GIGHA HERTZ   eV ou cm-1
% Gh=retdb(ld,'ld_2_Gh');(ou Gh=retdb(ld);)
% THz=retdb(ld,'ld_2_Thz');
% eV=retdb(ld,'ld_2_eV');
% cmm1=retdb(ld,'ld_2_cmm1');
% rad_s=retdb(ld,'ld_2_rad_s');
% ld=retdb(eV,'eV_2_ld');
% ld=retdb(Gh,'Gh_2_ld');
% ld=retdb(THz,'THz_2_ld');
% ld=retdb(cmm1,'cmm1_2_ld');
% ld=retdb(rad_s,'ld_2_rad_s');
% db=retdb(attenuation_en_energie,'attenuation_en_energie_2_db');
% attenuation_en energie=retdb(db,'db_2_attenuation_en energie');
%
% Definition des db  db=10*log10(Energie2/Energie1)
% See also: RETCONSTANTES



if nargin<2;db=2.99792458e5./n;return;end;
if isnumeric(ld);db=imag(n).*(40000*pi./(log(10)*ld));return;end
ct=retconstantes;

q=ct.e;h=ct.h;6.6261E-34;c=ct.c;
switch ld;
case {'ld_2_THz','THz_2_ld'};db=(1.e-6*c)./n;
case {'ld_2_Gh','Gh_2_ld'};db=(1.e-3*c)./n;
case {'ld_2_eV','eV_2_ld'};db=(h*c/(1.e-6*q))./n;
case {'ld_2_cmm1','cmm1_2_ld'};db=1.e4./n;
case {'ld_2_rad_s','rad_s_2_ld'};db=2*pi*c./(1.e-6*n);
case 'Gh_2_eV';db=(h/(1.e-9*q)).*n;
case 'eV_2_Gh';db=(1.e-9*q/h).*n;
case 'attenuation_en_energie_2_db';db=10*log10(n);
case 'db_2_attenuation_en energie';db=10.^(n/10);
%otherwise;db=imag(n).*(40000*pi./(log(10)*ld));
end;