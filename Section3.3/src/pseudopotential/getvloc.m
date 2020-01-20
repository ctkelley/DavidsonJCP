function [vion,rho] = getvloc(mol)
% GETVLOC returns the local pseudo potential of the molecule.
%    [vion,rho] = GETVLOC(mol) calculate the local pseudo potential of the
%    molecule based on the local pseudo potential of each atom in the
%    molecule in G-space. And return the corresponding atomic charge
%    density as well.
%
%    See also vnl2g, vloc2g, getvnl. 

%  Copyright (c) 2015-2016 Yingzhou Li and Chao Yang,
%                          Stanford University and Lawrence Berkeley
%                          National Laboratory
%  This file is distributed under the terms of the MIT License.

n1       = mol.n1;
n2       = mol.n2;
n3       = mol.n3;
vol      = mol.vol;
atoms    = mol.atoms;

vion  = zeros(n1,n2,n3);
rho   = zeros(n1,n2,n3);

if isempty(atoms)
    return;
end

xyzlist  = mol.xyzlist;
ntypes = length(mol.natoms);
alist = mol.alist;

vloc3dg = zeros(n1,n2,n3);
rho3dg  = zeros(n1,n2,n3);

ppvar = mol.ppvar;
vlocg = ppvar.vlocg;
rhog  = ppvar.rhog;

ecut2  = mol.ecut2;
grid2 = Ggrid(mol,ecut2);
inz    = grid2.idxnz;
gkx2   = grid2.gkx;
gky2   = grid2.gky;
gkz2   = grid2.gkz;

% vectorize the original for i = 1:nnz loop
for itype = 1:ntypes
    index  = alist == itype;
    xyzmat = xyzlist(index,:)';
    phmat  = [gkx2 gky2 gkz2]*xyzmat;
    %
    % ccvec is the structure factor used to account for
    % the translation of the atomic position
    %
    ccvec  = sum(exp(-1i.*phmat),2);
    vloc3dg(inz)  = vloc3dg(inz) + vlocg(:,itype).*ccvec;
    rho3dg(inz)  = rho3dg(inz) + rhog(:,itype).*ccvec;
end

% vloc3dg = conj(vloc3dg);
% rho3dg = conj(rho3dg);

%
% get rho in real space 
vion  = real(ifft3(vloc3dg)*n1*n2*n3)*4/vol*pi;
rho   = real(ifft3(rho3dg)*n1*n2*n3)/vol;

end
