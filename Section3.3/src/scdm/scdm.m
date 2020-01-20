function [Phi,Q] = scdm(Psi,varargin)
% SCDM Selected columns of the density matrix.
%    [Phi,Q] = SCDM(Psi) conducts compressed representation of Kohn?Sham
%    orbitals via selected columns of the density matrix. Phi returns
%    the sparse representation of the wavefunction and Q is the
%    corresponding rotation matrix.
%
%    [Phi,Q] = SCDM(Psi,methodstr) conducts compressed representation of
%    Kohn?Sham orbitals via selected columns of the density matrix with
%    the given method. methodstr can be either 'plain', 'rand' or 'fast'
%    which correspond to two different techniques. The default tolerance
%    is 5e-2.
%
%    [Phi,Q] = SCDM(Psi,methodstr,tol) conducts compressed
%    representation of Kohn?Sham orbitals via selected columns of the
%    density matrix with the given method and tolerance.
%
%   See also scdm_fast, scdm_rand.

%  Copyright (c) 2016-2018 KSSOLV Team
%  This file is distributed under the terms of the MIT License.

% set parameters
if nargin < 3
    tol = 5e-2;
else
    tol = varargin{2};
end

if nargin < 2
    methodstr = 'plain';
else
    methodstr = varargin{1};
end

switch lower(methodstr)
    case 'plain'
        [Q,~,~] = qr(Psi',0);
        Phi = Psi*Q;
    case 'rand'
        [Phi,Q] = scdm_rand(Psi);
    case 'fast'
        [Phi,Q] = scdm_fast(Psi,tol);
    otherwise
        error('Invalid methodstr input');
end

end