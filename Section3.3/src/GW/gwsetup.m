function ksinfo = gwsetup(mol)
% usage: set up a bunch of things for GW calculation
% ksinfo = gwsetup(mol)
% ksinfo is a structure that contains the following entries:
%
% nv    number of occupied states (exclude spin)
% Z     contains all eigenvectors of the KS Hamiltonian
% ev    contains the corresponding KS eigenvalues
% F     is the Fourier transform object (non-square) that maps planewave coefficients 
%       to Kohn-Sham orbitals in reals space
% coulG is the Coulomb potential in Fourier space (diagonal of the Coulomb matrix)
% vol   volume of the supercell
% ntot  Total number of grid points in 3D
%
%

grid = Ggrid(mol);
gkk = grid.gkk; % get a vector |G|^2 within the Ecut limit
idxnz = grid.idxnz; % get the position of gkk within a cube

% construct the discrete Fourier transformation matrix which n123 by ng
ksinfo.F = KSFFT(mol);
FF = ksinfo.F;
[ng,nr]=size(ksinfo.F);
if (ng ~= length(idxnz)) 
   fprintf('ng = %d, length(idxnz) = %d\n', ng, length(idxnz));
   error('src/GW/gwsetup.m:inconsistent ng and idxnz size');
end;

% run SCF to get ground state information
%
kssolvpptype('pz-hgh','UPF');
ksinfo.nv = get(mol,'nel')/2; % number of occupied states
[mol, H, X, info] = scf(mol);
Hmat = ham2mat(H); % turn it into a matrix
Hmat = (Hmat+Hmat')/2; % symmetrize to keep eigenvalues real 
% (causes problem for inversion symmetry?)
[ksinfo.Z,D] = eig(Hmat);
[ksinfo.ev,id]=sort(real(diag(D)));
%ksinfo.ev = ksinfo.ev*2; % convert to Rydberg
ksinfo.Z = ksinfo.Z(:,id);
ksinfo.ntot =  get(mol,'n1')*get(mol,'n2')*get(mol,'n3');
ksinfo.vol = det(get(mol,'supercell'));
%

% construct the bare Coulomb in reciprocal space
coulG = zeros(ng,1);
scell = mol.supercell;
amin = min([norm(scell(:,1)) norm(scell(:,2)) norm(scell(:,3))])/2; % not quite right
for j = 1:length(idxnz)
   if ( abs(gkk(j)) ~= 0 )
      coulG(j) = 4.0*pi/gkk(j); 
      coulG(j) = coulG(j)*(1-cos(sqrt(gkk(j))*amin));
   else
      %coulG(j) = 4.0*pi*amin^2/2;
      coulG(j) = 0.0;
   end;
end;
ksinfo.coulG = coulG;
ksinfo.coulG0 = 4.0*pi*amin^2/2;
