function [cry,BX,info] = nscf4c(cry,options)
% NSCF4C Non Self Consistent Field iteration for crystal.
%    [cry,BX,info] = NSCF4C(cry,options) find the ground state minimum
%    total energy and the corresponding wave functions for cry. The initial
%    density must be provided in options.
%
%   See also scf4m, scf4c, dcm, trdcm.

%  Copyright (c) 2016-2017 Yingzhou Li and Chao Yang,
%                          Stanford University and Lawrence Berkeley
%                          National Laboratory
%  This file is distributed under the terms of the MIT License.

if nargin < 2 || isempty(options.rho0)
    error('The initial density is not provided');
end

% Set timer
tstart  = cputime;

% Initialize input variables
rho        = options.rho0;
BX0        = options.X0;
nkpts      = cry.nkpts;
temperature= cry.temperature;
Tbeta      = temperature*8.6173324e-5/13.6;

% Initialize Hamiltonian, Wavefun, and Preconditioners
[cry,BH,BX,BHprec,noccs] = iterinit(cry,rho,BX0);

nBXcols    = ncols(BX);

% Initialize output variables
ev         = zeros(sumel(nBXcols),1);

fprintf('Beging NSCF4C calculation for %s...\n',cry.name);

idx = 0;
for ik = 1:nkpts
    idx = idx(end) + (1:nBXcols(ik));
    [X, ev(idx)] = updateX(cry, BH{ik}, BX{ik}, BHprec{ik}, options);
    BX{ik} = X;
end

[occs,~] = getocc(ev,noccs,Tbeta);

BX = assignoccs(BX,occs);

info.Eigvals = ev;

timetot = cputime - tstart;
fprintf('Total time used = %20.3e\n', timetot);

end