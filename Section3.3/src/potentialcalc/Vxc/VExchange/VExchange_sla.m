function [vx,ux] = VExchange_sla(rs)
% VEXCHANGE_SLA Slater exchange with alpha=2/3.
%    [vx,ux] = VEXCHANGE_SLA(rs) returns the Slater exchange with alpha=2/3
%    of the rs = (3/4*pi/rho)^(1/3).
%
%   See also exRef.

%  Copyright (c) 2016-2017 Yingzhou Li and Chao Yang,
%                          Stanford University and Lawrence Berkeley
%                          National Laboratory
%  This file is distributed under the terms of the MIT License.

% alpha = 2/3;
% f = -9/8*(3/2/pi)^(2/3);
% falpha = f*alpha;
falpha = -0.458165293283143;

vx = 4/3*falpha./rs;
ux = falpha./rs;

end