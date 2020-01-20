function [d, VA, WA, VB, WB, Hbse] = ...
    getkernel(ksinfo, omega, eta, nvbands, ncbands, W)
%
%  [d, VA, WA, VB, WB] = ...
%      construct_kernel(ksinfo, omega, eta, nbands, W);
%
% Construct BSE kernel.
%
% ksinfo is a structure that contains ground state calculation results
% and parameters of the molecules
%   nv      --- number of valence states
%   Z       --- contains the eigenvecgtors from the KSDFT calculation
%   ev      --- contains the corresponding eigenvalues
%   vol     --- volume of the unit cell
%   ntot    --- total number of grid points on which the wavefunction is sampled
% 
% omega   --- frequency at which the BSE kernel is being evaluated
% eta     --- Lorentzian broadening factor that turns a Dirac delta into a 
%             a smoother peak
% nvbands --- the number of valence bands (KS orbitals) from the Fermi 
%             level included in the kernel calculation. 
% ncbands --- the number of empty bands (KS orbitals) included in 
%             the kernel calculation
%
% d       --- a vector of length nvbands*ncbands that contains the
%             energy difference between valence and empty bands.
%             To be specific, for an (iv,ic) pair
%              d((iv-1)*ncbands+ic) = ev(nv+ic) - ev(nv-nvbands+iv);
% VA      --- an nc x nv x nc x nv tensor that contains the exchange
%             integrals  
%                VA(i,j,ip,jp) = < phi(i)*phi(j) | v | phi(ip)*phi(jp) >
% VB      --- similar to VA, but with the expression
%                VB(i,j,ip,jp) = < phi(i)*phi(j) | v | conj(phi(ip)*phi(jp)) >
%
% WA      --- a nc x nc x nv x nv tensor that contains the direct integrals
%                WA(i,j,ip,jp) = < phi(i)*phi(j) | W | phi(ip)*phi(jp) >
%             where i, j are valence bands and ip, jp are empty bands
%
% WB      --- a nv x nc x nv x nc tensor that contains the direct integrals
%                WB(i,j,ip,jp) = < phi(i)*phi(j) | W | phi(ip)*phi(jp) >
%             where i, ip are valence bands and j, jp are empty bands
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
Dcoul = spdiags(ksinfo.coulG,0,ng,ng);

% Normalize the wavefunction in the Fourier space
% 1/Vol sum_G |Psi(G)|^2 = 1
for iv = 1 : ng
  Z(:,iv) = Z(:,iv) * sqrt(vol) / (norm(Z(:,iv)));
end

% Ordering:
% [v1*c1, ..., v1*cn, v2*c1, ..., v2*cn, ..., vk*c1, ..., vk*cn]
%d = zeros(nc*nv, 1);
%for iv = 1:nv,
%  d((iv-1)*nc + (1:nc)) = ev(nv+1:n) - ev(iv);
%end
d = zeros(ncbands*nvbands, 1);
for iv = 1:nvbands
  for ic = 1:ncbands
    d((iv-1)*ncbands+ic) = ev(nv+ic) - ev(nv-nvbands+iv);
  end
end

if (ng*nvbands*ncbands*16 + nr*16*2 < 1e9) 
   % use BLAS3 if there is enough memory
   % Ordering:
   % [v1*c1, ..., v1*cn, v2*c1, ..., v2*cn, ..., vk*c1, ..., vk*cn]
   %psivc = [];
   psivc = zeros(ng, ncbands*nvbands);
   for iv = 1:nvbands
     psivr = conj(F'*Z(:,nv-nvbands+iv));
     for ic = 1:ncbands
       psivc(:,(iv-1)*ncbands+ic) = F*((F'*Z(:,nv+ic)).*psivr);
     end
   end
   % Ordering:
   % [v1*v1, ..., v1*vk, v2*v1, ..., v2*vk, ..., vk*v1, ..., vk*vk]
   psivv = zeros(ng, nvbands*nvbands);
   for iv = 1:nvbands
     psivr = conj(F'*Z(:,nv-nvbands+iv));
     for jv = 1:nvbands
       psivv(:,(iv-1)*nvbands+jv) = F*((F'*Z(:,nv-nvbands+jv)).*psivr);
     end
   end
   % Ordering:
   % [c1*c1, ..., c1*cn, c2*c1, ..., c2*cn, ..., cn*c1, ..., cn*cn]
   psicc = zeros(ng, ncbands*ncbands);
   for ic = 1:ncbands
     psicr = conj(F'*Z(:,nv+ic));
     for jc = 1:ncbands
       psicc(:,(ic-1)*ncbands+jc) = F*((F'*Z(:,nv+jc)).*psicr);
     end
   end
  

%   psivcR = F'*psivc;
%   psivvR = F'*psivv;
%   psiccR = F'*psicc;

%size(psivcR)
%[UUvc, SSvc, VVvc] = svd(psivcR, 0);

%size(psivvR)
%[UUvv, SSvv, VVvv] = svd(psivvR, 0);

%size(psiccR)
%[UUcc, SScc, VVcc] = svd(psiccR, 0);



   
   VA = psivc'*(Dcoul*psivc)/(4.0*pi);
   VB = psivc'*(Dcoul*conj(psivc))/(4.0*pi);
   
   WA0 = psicc'*W*psivv/(4.0*pi);
   WB0 = psivc'*W*psivc/(4.0*pi);
   
   WA0head = psicc(1,:)'*W(1,1)*psivv(1,:)/(4.0*pi);
   WB0head = psivc(1,:)'*W(1,1)*psivc(1,:)/(4.0*pi);
   
   WA0wing1 = psicc(1,:)'*W(1,2:end)*psivv(2:end,:)/(4.0*pi);
   WB0wing1 = psivc(1,:)'*W(1,2:end)*psivc(2:end,:)/(4.0*pi);
   
   WA0wing2 = psicc(2:end,:)'*W(2:end,1)*psivv(1,:)/(4.0*pi);
   WB0wing2 = psivc(2:end,:)'*W(2:end,1)*psivc(1,:)/(4.0*pi);
   
   WA0body = psicc(2:end,:)'*W(2:end,2:end)*psivv(2:end,:)/(4.0*pi);
   WB0body = psivc(2:end,:)'*W(2:end,2:end)*psivc(2:end,:)/(4.0*pi);

else
   % use a lot less memory, but slow
   for ivp = 1:nvbands
      psivpr = conj(F'*Z(:,nv-nvbands+ivp));
      for icp = 1:ncbands
         psivpcp = F*((F'*Z(:,nv+icp)).*psivpr);
         for iv = 1:nvbands
            psivr = conj(F'*Z(:,nv-nvbands+iv));
            for ic = 1:ncbands;
               psivc = F*((F'*Z(:,nv+ic)).*psivr);
               VA(ic,iv,icp,ivp) = psivc'*(Dcoul*psivpcp)/(4.0*pi);
               VB(ic,iv,icp,ivp) = psivc'*(Dcoul*conj(psivpcp))/(4.0*pi);
               WB0(ic,iv,icp,ivp) = psivc'*(W*psivpcp)/(4.0*pi);
            end
         end
      end
   end

   for ivp = 1:nvbands
      psivrp = conj(F'*Z(:,nv-nvbands+ivp));
      for iv = 1:nvbands
         psivv = F*((F'*Z(:,nv-nvbands+iv)).*psivrp);
         for icp = 1:ncbands
            psicrp = conj(F'*Z(:,nv+icp));
            for ic = 1:ncbands
               psicc = F*((F'*Z(:,nv+ic)).*psicrp);
               WA0(icp,ic,ivp,iv) = psicc'*(W*psivv)/(4.0*pi);
            end
         end
      end
   end
end

for iv = 1:nvbands,
  for jc = 1:ncbands,
    for ivp = 1:nvbands,
      for jcp = 1:ncbands,
        rl = (iv-1)*ncbands+jc;
        cl = (ivp-1)*ncbands+jcp;
        rr = (jc-1)*ncbands+jcp;
        cr = (iv-1)*nvbands+ivp;
        WA(rl, cl) = WA0(rr, cr);
        WAhead(rl, cl) = WA0head(rr, cr);
        WAwing1(rl, cl) = WA0wing1(rr, cr);
        WAwing2(rl, cl) = WA0wing2(rr, cr);
        WAbody(rl, cl) = WA0body(rr, cr);
        rl = (iv-1)*ncbands+jc;
        cl = (ivp-1)*ncbands+jcp;
        rr = (ivp-1)*ncbands+jc;
        cr = (iv-1)*ncbands+jcp;
        WB(rl, cl) = WB0(rr, cr);
        WBhead(rl, cl) = WB0head(rr, cr);
        WBwing1(rl, cl) = WB0wing1(rr, cr);
        WBwing2(rl, cl) = WB0wing2(rr, cr);
        WBbody(rl, cl) = WB0body(rr, cr);
      end
    end
  end
end

%inveps0 = inveps(1);
%w_eff = vcoul0 * inveps0/(8*pi);
bsemat_fac = -8.00*pi/vol; 

%WA = WAhead + WAwing1 + WAwing2 + WAbody; 
%WB = WBhead + WBwing1 + WBwing2 + WBbody; 

Kd = bsemat_fac * WA;
Kx = (-2.00 * bsemat_fac) * VA;

Hbse = Kd + Kx;

KK = zeros(ncbands*nvbands, ncbands*nvbands);
VB = KK;
WB = KK;
