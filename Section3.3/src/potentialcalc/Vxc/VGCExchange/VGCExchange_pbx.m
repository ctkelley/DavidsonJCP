function [v1gcx,v2gcx,ugcx] = VGCExchange_pbx(rho,grho2)
% VGCEXCHANGE_PBX PBE exchange gradient correction.
%    [v1gcx,v2gcx,ugcx] = VGCEXCHANGE_PBX(rho,grho2) returns the PBE
%    exchange gradient correction of the rho and the gradient square of
%    rho.
%
%   See also exRef.

%  Copyright (c) 2016-2017 Yingzhou Li and Chao Yang,
%                          Stanford University and Lawrence Berkeley
%                          National Laboratory
%  This file is distributed under the terms of the MIT License.

c1 = 0.75/pi;
c2 = 3.093667726280136;
c5 = 4/3;

k = 0.804;
mu = 0.21951;

agrho = sqrt(grho2);
kf = c2*rho.^(1/3);
dsg = 0.5./kf;
s1 = agrho.*dsg./rho;

% Energy
f2 = 1 + s1.*s1*mu/k;
fx = k - k./f2;

exunif = -c1*kf;
ugcx = exunif.*fx;

% Potential
dxunif = exunif/3;

dfx = 2 * mu * s1 ./ f2.^2;
v1gcx = ugcx + dxunif.*fx - exunif.*dfx*c5.*s1;
v2gcx = exunif.*dfx.*dsg./agrho;

end