function ebands = eband(cry,rho,nbnd,stpt,edpt)
% EBAND calculates the energy band.
%   ebands = EBAND(cry,rho,nbnd,kpts) calculates the smallest nbnd
%   eigenvalues of a Hamiltonian reconstructed by cry and rho at the given
%   k-points. ebands is a nkpts by nbnd matrix and each row of ebands is
%   the eigenvalues of the corresponding k-points.
%
%   ebands = EBAND(cry,rho,nbnd,stpt,edpt) calculates the smallest nbnd
%   eigenvalues of a Hamiltonian reconstructed by cry and rho at the given
%   line of k-points with stpt and edpt being the starting and ending
%   points.
%
%   See also plotband.

%  Copyright (c) 2016-2017 Yingzhou Li, Lin Lin and Chao Yang,
%                          Stanford University,
%                          University of California, Berkeley
%                          and Lawrence Berkeley National Laboratory
%  This file is distributed under the terms of the MIT License.

if nargin < 5
    kpts = stpt;
else
    gap = 0.001;
    len = sqrt(sum((edpt-stpt).^2));
    kpts = (0:gap:len)'*(edpt-stpt);
    kpts = kpts + repmat(stpt,size(kpts,1),1);
end

crykpts = set(cry,'kpts',kpts);
crykpts = finalize(crykpts);

opt = setksopt();
opt.rho0 = rho;
opt.X0 = genX0(crykpts,nbnd*ones(size(kpts,1),1));

[~,~,infokpts] = nscf4c(crykpts,opt);

ebands = reshape(infokpts.Eigvals, nbnd, [])';

end