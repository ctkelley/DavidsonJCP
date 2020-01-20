function cry = set(cry,varargin)
% CRYSTAL/SET Set function for crystal class
%    cry = SET(cry,str1,field1,str2,field2,...) returns a crystal class of
%    the given fields with respect to the name strings.
%
%    See also Molecule, Crystal.

%  Copyright (c) 2016-2017 Yingzhou Li and Chao Yang,
%                          Stanford University and Lawrence Berkeley
%                          National Laboratory
%  This file is distributed under the terms of the MIT License.

nvar = length(varargin);
if mod(nvar,2) == 1
    error('Wrong input for Crystal.set');
end

for it = 1:2:nvar
    attr_name = varargin{it};
    value     = varargin{it+1};
    
    if strcmpi(attr_name,'atomlist')
        [~,IA,cry.alist] = unique([value.anum]);
        cry.atoms = value(IA);
        continue;
    end
    
    if strcmpi(attr_name,'autokpts')
        nkx = value(1);
        nky = value(2);
        nkz = value(3);
        if numel(value) ==3
            skx = 0;
            sky = 0;
            skz = 0;
        else
            skx = value(4);
            sky = value(5);
            skz = value(6);
        end
        [I,J,K] = ndgrid((0:nkx-1),(0:nky-1),(0:nkz-1));
        pregkx = (I(:)+skx)/nkx;
        pregky = (J(:)+sky)/nky;
        pregkz = (K(:)+skz)/nkz;
        cry.kpts = [pregkx pregky pregkz];
        continue;
    end
    
    cry.(attr_name) = value;
end

end
