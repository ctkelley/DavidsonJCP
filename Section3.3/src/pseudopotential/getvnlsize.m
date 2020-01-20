function [siz,offset] = getvnlsize(mol)
% GETVNLSIZE returns the size of non-local pseudo potential.
%    [vnlmat,vnlsign] = GETVNLSIZE(mol) returns the size of the non-local
%    pseudo potential for each atom of the molecule.
%
%    See also getvloc, getvnl.

%  Copyright (c) 2015-2016 Yingzhou Li and Chao Yang,
%                          Stanford University and Lawrence Berkeley
%                          National Laboratory
%  This file is distributed under the terms of the MIT License.


ppvar = mol.ppvar;
alist = mol.alist;
ntypes = length(mol.atoms);
na = length(alist);

siz = zeros(na,1);
for itype = 1:ntypes
    totall = 0;
    for it = 1:length(ppvar.lll{itype})
        ll = ppvar.lll{itype}(it);
        totall = totall + 2*ll+1;
    end
    index  = find(alist == itype);
    for ita = index'
        siz(ita) = totall;
    end
end

offset = cumsum(siz);
offset = [0; offset(1:end-1)];

end