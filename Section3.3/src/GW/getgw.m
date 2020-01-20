function [Esx, Ech, Ex] = getgw(ksinfo, eta, nvbands, ncbands, chi0)
%
%
% ksinfo is a structure that contains ground state calculation results
% and parameters of the molecules
%   nv      --- number of valence states
%   Z       --- contains the eigenvecgtors from the KSDFT calculation
%   ev      --- contains the corresponding eigenvalues
%   vol     --- volume of the unit cell
%   ntot    --- total number of grid points on which the wavefunction is sampled
% 
% eta     --- Lorentzian broadening factor that turns a Dirac delta into a 
%             a smoother peak
% nvbands --- the number of valence bands (KS orbitals) from the Fermi 
%             level included in the kernel calculation. 
% ncbands --- the number of empty bands (KS orbitals) included in 
%             the kernel calculation
%
%
%
nv   = ksinfo.nv;
Z    = ksinfo.Z;
ev   = ksinfo.ev;
F    = ksinfo.F;
vol  = ksinfo.vol;
ntot = ksinfo.ntot;
%
% 1:nv = occupied, nv+1:n = unoocupied
im = sqrt(-1);
[ng,nr]=size(F);

% Normalize the wavefunction in the Fourier space
% 1/Vol sum_G |Psi(G)|^2 = 1
for iv = 1 : ng
  Z(:,iv) = Z(:,iv) * sqrt(vol) / (norm(Z(:,iv)));
end

ng = size(chi0,1);
I = eye(ng);
vol  = ksinfo.vol;

Dcoul0 = spdiags(ksinfo.coulG,0,ng,ng);

%ksinfo.coulG(1) = ksinfo.coulG0;
epsilon = geteps(chi0, ksinfo.coulG);

inveps = inv(epsilon);
nbands = nvbands + ncbands;

ksinfo.coulG(1) = ksinfo.coulG0;
Dcoul   = spdiags(ksinfo.coulG,0,ng,ng);
W1 = epsilon\Dcoul;
W2 = W1 - Dcoul;

Esx = zeros(nbands, nbands);
Ech = zeros(nbands, nbands);
Ex = zeros(nbands, nbands);

for n1 = 1:nbands
  psin1r = conj(F'*Z(:,n1));
  for n2 = 1:nbands
    psin2r = conj(F'*Z(:,n2));
    for iocc = 1:nv
      psioccr = F'*Z(:,iocc);
      psioccn1 = F*(psioccr.*psin1r);
      psioccn2 = F*(psioccr.*psin2r);
      %Esx(n1, n2) = Esx(n1, n2) + psioccn1'*inveps*Dcoul*psioccn2;
      Esx(n1, n2) = Esx(n1, n2) + psioccn1'*(W1*psioccn2);
      Ech(n1, n2) = Ech(n1, n2) + psioccn1'*(W2*psioccn2);
      Ex(n1, n2) = Ex(n1, n2) + psioccn1'*(Dcoul*psioccn2);
    end  
  end
end

%fac = 8.00*pi/vol 
fac = 1.00/vol 

Esx = 0.0 - fac * Esx;
Ech = 0.5 * fac * Ech;
Ex = 0.0 - fac * Ex;
