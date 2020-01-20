function entropy = getEntropy(occ,Tbeta)
% GETENTROPY calculates entropy of the statues.
%    entropy = GETENTROPY(occ,Tbeta) calculates the entropy of the
%    occupation status and bete temperature. The calculation is very
%    standard, i.e., entropy = sum_i occ(i)log occ(i).
%
%   See also scf4m.

%  Copyright (c) 2016-2017 Yingzhou Li and Chao Yang,
%                          Stanford University and Lawrence Berkeley
%                          National Laboratory
%  This file is distributed under the terms of the MIT License.

if Tbeta < eps
    entropy = 0;
else
    occ = occ.*(occ>eps) + (1-occ).*(occ<eps);
    entropy = sum(occ.*log(occ)) * Tbeta;
end

end