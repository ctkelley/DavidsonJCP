function [v1gcc,v2gcc,ugcc] = VGCCorrelation_pbc(rho,grho2,ec,vc)
% VGCCORRELATION_PBC PBE correlation gradient correction.
%    [v1gcc,v2gcc,ugcc] = VGCCORRELATION_PBC(rho,grho2) returns the PBE
%    correlation gradient correction of the rho and the gradient square of
%    rho.
%
%   See also exRef.

%  Copyright (c) 2016-2017 Yingzhou Li and Chao Yang,
%                          Stanford University and Lawrence Berkeley
%                          National Laboratory
%  This file is distributed under the terms of the MIT License.

ga = 0.031091;
be = 0.066725;
xkf = 1.919158292677513;
xks = 1.128379167095513;

if nargin < 3
    rs = ((3/(4*pi))./rho).^(1/3);
    [vc,ec] = VCorrelation_pw(rs);
end

ks = xks.*sqrt(xkf./rs);
t2 = grho2./(2*ks.*rho).^2;
expe = exp(-ec/ga);
af = be/ga*(1./(expe-1));
bf = expe.*(vc-ec);
y = af.*t2;
xy = (1+y)./(1+y+y.*y);
qy = y.*y.*(2+y)./(1+y+y.*y).^2;
s1 = 1+be/ga*t2.*xy;
h0 = ga*log(s1);
dh0 = be*t2./s1.*(-7/3*xy-qy.*(af.*bf/be-7/3));
ddh0 = be./(2*ks.*ks.*rho).*(xy-qy)./s1;

ugcc = h0;
v1gcc = h0 + dh0;
v2gcc = ddh0;

end