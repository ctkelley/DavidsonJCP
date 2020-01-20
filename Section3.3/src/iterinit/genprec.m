function prec = genprec(H)
% GENPREC generates the preconditioner.
%    X0 = GENPREC(H) generates the preconditioner for the Hamiltonian,
%    which is in the format of Wavefun.
%
%    BX0 = GENPREC(BH) generates the preconditioners for the Bloch
%    Hamiltonian, which is in the format of BlochWavefun.
%
%   See also Ham, BlochHam, Wavefun, BlochWavefun.

%  Copyright (c) 2016-2017 Yingzhou Li and Chao Yang,
%                          Stanford University and Lawrence Berkeley
%                          National Laboratory
%  This file is distributed under the terms of the MIT License.

idxnz = H.idxnz;
n1 = H.n1;
n2 = H.n2;
n3 = H.n3;
if isa( H, 'BlochHam' )
    nkpts = H.nkpts;
    gkincell = H.gkincell;
else
    nkpts = 1;
    gkincell = {H.gkin};
end

p = cell(nkpts,1);
for ik = 1:nkpts
    X  = gkincell{ik};
    Y  = 27.0 + X.*(18.0 + X.*(12.0 + 8.0*X));
    p{ik}  = Y./(Y + 16.0*X.^4);
end;

if isa( H, 'BlochHam' )
    prec = BlochWavefun(p,n1,n2,n3,idxnz,H.wks);
else
    prec = Wavefun(p{1},n1,n2,n3,idxnz);
end

end