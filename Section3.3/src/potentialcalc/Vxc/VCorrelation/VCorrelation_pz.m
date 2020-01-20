function [vc,uc] = VCorrelation_pz(rs)
% VCORRELATION_PZ Perdew-Zunger correlation.
%    [vc,uc] = VCORRELATION_PZ(rs) returns the Perdew-Zunger correlation of
%    the rs = (3/4*pi/rho)^(1/3).
%
%   See also exRef.

%  Copyright (c) 2016-2017 Yingzhou Li and Chao Yang,
%                          Stanford University and Lawrence Berkeley
%                          National Laboratory
%  This file is distributed under the terms of the MIT License.

a = 0.0311;
b = -0.048;
c = 0.0020;
d = -0.0116;
gc = -0.1423;
b1 = 1.0529;
b2 = 0.3334;

vc = zeros(size(rs));
uc = zeros(size(rs));

% high density formula
idxl = rs < 1;
rsl = rs(idxl);
lnrs = log (rsl);
uc(idxl) = a*lnrs + b + c*rsl.*lnrs + d*rsl;
vc(idxl) = a*lnrs + (b-a/3) + 2/3*c*rsl.*lnrs + (2*d-c)/3*rsl;

% interpolation formula
idxg = rs >= 1;
rsg = rs(idxg);
rs12 = sqrt(rsg);
ox = 1 + b1*rs12 + b2*rsg;
dox = 1 + 7/6*b1*rs12 + 4/3*b2*rsg;
uc(idxg) = gc./ox;
vc(idxg) = uc(idxg).*dox./ox;

end
