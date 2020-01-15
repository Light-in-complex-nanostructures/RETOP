function [v,vv,v6]=retversion(v6);
% function [v,vv,v6]=retversion(v6);
% version pour imposer le format d' ecriture de retio, retversion('v7'), retversion('v7.3'),
% si v6(chaine de caracteres)  est en entree, en sortie v6=['-',v6]
% vv est un nombre croissant avec les versions
% 7.3.0.267  ==> v=7  vv=7.003
% 7.7.0.471  ==> v=7  vv=7.007
% 7.13       ==> v=7  vv=7.013
% 8.2        ==> v=8  vv=8.002
persistent vv6 VV V;
if isempty(V);
vv=version;v=sscanf(vv,'%d');%vv=sscanf(vv,'%f');vv=vv(1)
points=findstr(vv,'.');vv=str2num(vv(1:points(1)))+.001*str2num(vv(points(1)+1:points(2)-1));V=v;VV=vv;
end;
v=V;vv=VV;
if isempty(vv6);
    if (v>6);
    if vv>=7.7;vv6='-v7';else;vv6='-v6';end;
    else;vv6='';
    end;
end;
if nargout>2;
if nargin<1;v6=vv6;	
else;
if vv>=7.3;v6=['-',v6];vv6=v6;end;
end;
end;	
	
	
