function rho = getcharge(mol,X,occ)
% GETCHARGE calcuates the charge density.
%    rho = GETCHARGE(mol,X) calculates the charge density generated by X.
%
%    rho = GETCHARGE(mol,X,occ) calculates the charge density generated by
%    X together with the occupation rate occ.
%
%    rho = GETCHARGE(cry,BX) calculates the charge density generated by
%    the Bloch wave function BX.
%
%    rho = GETCHARGE(cry,BX,occs) calculates the charge density generated
%    by the Bloch wave function BX together with the occupation rate occs.
%
%   See also scf, dcm, trdcm.

%  Copyright (c) 2016-2017 Yingzhou Li and Chao Yang,
%                          Stanford University and Lawrence Berkeley
%                          National Laboratory
%  This file is distributed under the terms of the MIT License.

if (nargin == 2)
    occ = ones(1,sumel(ncols(X)));
elseif iscell(occ)
    occ = cell2mat(occ);
else
    occ = occ(:)';
end

n1    = mol.n1;
n2    = mol.n2;
n3    = mol.n3;
n123  = n1*n2*n3;
nspin = mol.nspin;

if isa(mol,'Crystal')
    rho = zeros(n123,1);
    nXcols = ncols(X);
    idx = 0;
    for ik = 1:mol.nkpts
        idx = idx(end)+(1:nXcols(ik));
        psir = ifft3(X{ik})*n123;
        rho = rho+mol.wks(ik)*sum(repmat(occ(idx),n123,1).*abs(psir).^2,2);
    end
else
    psir = ifft3(X)*n123;
    rho = sum(repmat(occ,n123,1).*abs(psir).^2,2);
end
rho = reshape((2/nspin)*rho/mol.vol,n1,n2,n3);

end