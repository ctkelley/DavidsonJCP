function X0 = genX0(mol,nXcols)
% GENX0 generates the initial wave functions.
%    X0 = GENX0(mol,ncols) generates the initial wave functions for
%    molecule mol with ncols columns.
%
%    BX0 = GENX0(cry,ncols) generates the initial Bloch wave functions for
%    crystal cry with ncols columns for each k-points.
%
%   See also Molecule, Crystal, Wavefun, BlochWavefun.

%  Copyright (c) 2016-2017 Yingzhou Li and Chao Yang,
%                          Stanford University and Lawrence Berkeley
%                          National Laboratory
%  This file is distributed under the terms of the MIT License.

if nargin == 1
    if isa( mol, 'Crystal' )
        nXcols = mol.nel/2*mol.nspin*ones(mol.nkpts,1);
    else
        nXcols = mol.nel/2*mol.nspin;
    end
end

nkpts = numel(nXcols);

n1   = mol.n1;
n2   = mol.n2;
n3   = mol.n3;

grid = Ggrid(mol);
idxnz = grid.idxnz;

Qcell = cell(nkpts,1);

for ik = 1:nkpts
    psir = randn(n1,n2,n3,nXcols(ik));
    psif = reshape(fft3(psir),n1*n2*n3,[]);
    psif = psif(idxnz,:);
    [Qcell{ik},~]=qr(psif,0);
end

if isa( mol, 'Crystal' )
    X0 = BlochWavefun(Qcell,n1,n2,n3,idxnz,mol.wks);
else
    X0 = Wavefun(Qcell{1},n1,n2,n3,idxnz);
end

end
