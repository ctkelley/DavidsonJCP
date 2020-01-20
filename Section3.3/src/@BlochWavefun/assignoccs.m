function BX = assignoccs(BX,occs)
% BLOCHWAVEFUN/ASSIGNOCCS Assign occupation rate
%    BX = assignoccs(BX,occs) assigns the occupation rate to each of the
%    Bloch wave functions.
%
%    See also BlochWavefun.

%  Copyright (c) 2016-2017 Yingzhou Li and Chao Yang,
%                          Stanford University and Lawrence Berkeley
%                          National Laboratory
%  This file is distributed under the terms of the MIT License.

if iscell(occs)
    for it = 1:BX.nkpts
        BX.wavefuncell{it}.occ = occs{it};
    end
else
    idx = 0;
    for it = 1:BX.nkpts
        idx = idx(end) + (1:ncols(BX.wavefuncell{it}));
        BX.wavefuncell{it}.occ = occs(idx);
    end
end

end