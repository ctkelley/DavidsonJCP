function mol = set(mol,varargin)
% MOLECULE/SET Set function for molecule class
%    mol = SET(mol,str1,field1,str2,field2,...) returns a molecule class of
%    the given fields with respect to the name strings.
%
%    See also Molecule.

%  Copyright (c) 2016-2017 Yingzhou Li and Chao Yang,
%                          Stanford University and Lawrence Berkeley
%                          National Laboratory
%  This file is distributed under the terms of the MIT License.

nvar = length(varargin);
if mod(nvar,2) == 1
    error('Wrong input for Molecule.set');
end

for it = 1:2:nvar
    attr_name = varargin{it};
    value     = varargin{it+1};
    
    if strcmpi(attr_name,'atomlist')
        [~,IA,mol.alist] = unique([value.anum]);
        mol.atoms = value(IA);
        continue;
    end
    
    mol.(attr_name) = value;
end

end
