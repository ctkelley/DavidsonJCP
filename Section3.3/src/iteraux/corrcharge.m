function rho = corrcharge(mol,rho)
% CORRCHARGE correct the charge density.
%    rho = CORRCHARGE(mol,rho) correct the charge density, i.e., correct
%    negative density, renormalise density.
%
%   See also getcharge.

%  Copyright (c) 2016-2017 Yingzhou Li and Chao Yang,
%                          Stanford University and Lawrence Berkeley
%                          National Laboratory
%  This file is distributed under the terms of the MIT License.

n1    = mol.n1;
n2    = mol.n2;
n3    = mol.n3;
n123  = n1*n2*n3;
nel   = mol.nel;
charge = sumel(rho)*mol.vol/n123;

if abs(charge-nel)/nel > 1e-7
    warning(['Renormalize starting charge ' num2str(charge) ...
        ' to ' num2str(nel)]);
    if charge > 1e-8
        rho = rho*(nel/charge);
    else
        rho = nel/mol.vol*ones(size(rho));
    end
else

end