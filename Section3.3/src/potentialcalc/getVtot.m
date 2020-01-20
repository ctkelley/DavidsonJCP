function vtot = getVtot(mol,vion,vext,vhart,vxc)
%
% Usage: vtot = getvtot(mol,vion,vext,vhart,vxc);
%
% Purpose:
%   Compute the total potential energy and keep it band limited
%
% Input: 
%     mol   a Molecular object
%    vion   a 3D array that contains the local ionic potential
%    vext   a 3D array that contains external potential. In most
%           case, this should be an empty array
%    vhart  a 3D array that contains the Hartree (Coulomb) potential.
%    vxc    a 3D array that contains the exchange-correlation potential
%
% Output:
%    vtot   a 3D array that contains the total (local) potential
%

[n1,n2,n3]=size(vion);
ecut2   = mol.ecut2;
grid2  = Ggrid(mol,ecut2);
idxnz2  = grid2.idxnz;
vtot = vion + vhart + vxc;
if (~isempty(vext))
  vtot = vtot + vext;
end

% keep vtot band limited (sinc interpolation)
gm2 = zeros(n1,n2,n3);
gm2(idxnz2) = 1;
iz = gm2==0;
vfft = fftn(vtot);
vfft(iz) = 0;
vtot = ifftn(vfft);
