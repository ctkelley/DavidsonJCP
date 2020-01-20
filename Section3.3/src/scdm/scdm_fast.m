function [Phi,Q] = scdm_fast(Psi,tol)
% SCDM_FAST Selected columns of the density matrix with fast selection.
%    [Phi,Q] = SCDM_FAST(Psi) conducts compressed representation of
%    Kohn?Sham orbitals via selected columns of the density matrix with
%    fast selection with tol being the tolerance of indicator. Phi
%    returns the sparse representation of the wavefunction and Q is the
%    corresponding rotation matrix.
%
%   See also scdm, scdm_rand.

%  Copyright (c) 2016-2018 KSSOLV Team
%  This file is distributed under the terms of the MIT License.

[Phi_rand,Q_rand] = scdm_rand(Psi);
if isa(Phi_rand,'Wavefun')
    Psimat = Phi_rand.psi;
else
    Psimat = Phi_rand;
end

[~,ncol] = size(Psimat);

Psimat2 = Psimat.^2;
Idc = double(sparse(Psimat2 > tol*max(Psimat2,[],1)));
D = Idc'*Idc;

selected= [];
for k = 1:ncol
    cidx = find(D(:,k));
    ridx = find(sum(Idc(:,cidx)));
    [~,~,piv_loc] = qr(Psimat(ridx,cidx)',0);
    selected = union(selected,ridx(piv_loc(1:length(cidx))));
end

[Q_fast,~,~] = qr(Psimat(selected,:)',0);
Phi = Phi_rand*Q_fast;
Q = Q_rand*Q_fast;

end