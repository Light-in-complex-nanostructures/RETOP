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

function champ=extract_pointcoordinate_on_box(varargin);
coordonnees=varargin;
save prv_coordonnees coordonnees;
champ=zeros(length(varargin{1}),3*length(varargin)-3);
