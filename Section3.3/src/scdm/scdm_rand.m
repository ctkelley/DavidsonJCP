function [Phi,Q] = scdm_rand(Psi)
% SCDM_RAND Selected columns of the density matrix with random sampling.
%    [Phi,Q] = SCDM_RAND(Psi) conducts compressed representation of
%    Kohn?Sham orbitals via selected columns of the density matrix with
%    random sampling. Phi returns the sparse representation of the
%    wavefunction and Q is the corresponding rotation matrix.
%
%   See also scdm, scdm_fast.

%  Copyright (c) 2016-2018 KSSOLV Team
%  This file is distributed under the terms of the MIT License.

if isa(Psi,'Wavefun')
    Psimat = Psi.psi;
else
    Psimat = Psi;
end

w = sum(Psimat.^2,2);
[nrow,ncol] = size(Psimat);

for k = 2:5
    I = randsamppd(nrow,k*ncol,w);
    if cond(Psimat(I,:)) > 1e-8
        break;
    end
end

[Q,~,~] = qr(Psimat(I,:)',0);
Phi = Psi*Q;

end