function [vnlmatcell,vnlsigncell] = getvnlcell(cry)
% GETVNLCELL calculates the non-local pseudo potential for each k-points.
%    [vnlmatcell,vnlsigncell] = GETVNLCELL(cry,pseudovar) calculates the
%    non-local pseudo potential for each k-points in crystal.
%
%    See also VLOC2G, GETVNL.

%  Copyright (c) 2015-2016 Yingzhou Li and Chao Yang,
%                          Stanford University and Lawrence Berkeley
%                          National Laboratory
%  This file is distributed under the terms of the MIT License.

nkpts = cry.nkpts;
vnlmatcell = cell(nkpts,1);
vnlsigncell = cell(nkpts,1);
for ik = 1:nkpts
    [vnlmatcell{ik},vnlsigncell{ik}] = ...
        getvnl(cry,cry.kpts(ik,:));
end

end