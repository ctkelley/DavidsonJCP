function [d, VA, WA, VB, WB, Hbse] = ...
    getkernelisdf(ksinfo, omega, eta, nvbands, ncbands, W, vcrank_ratio, vvrank_ratio, ccrank_ratio)
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


%ISDF

% size of psivr: nr * nvbands; size of psicr: nr * ncbands
psivr = F'*Z(:,nv-nvbands+1:nv);
psicr = F'*Z(:,nv+1:nv+ncbands);

% size of vcrzeta_mu: nr * vcrank_mu
vcrank_mu = ceil(nvbands*ncbands*vcrank_ratio);
[vcrzeta_mu, vcind_mu] = isdf(psicr, conj(psivr), vcrank_mu);
vczeta_mu = F*vcrzeta_mu;
psivc_mu = prod_states(psicr(vcind_mu, :), conj(psivr(vcind_mu, :)));

% size of vvrzeta_mu: nr * vvrank_mu
vvrank_mu = ceil(nvbands*nvbands*vvrank_ratio);
[vvrzeta_mu, vvind_mu] = isdf(psivr, conj(psivr), vvrank_mu);
vvzeta_mu = F*vvrzeta_mu;
psivv_mu = prod_states(psivr(vvind_mu, :), conj(psivr(vvind_mu, :)));

% size of ccrzeta_mu: nr * ccrank_mu
ccrank_mu = ceil(ncbands*ncbands*ccrank_ratio);
[ccrzeta_mu, ccind_mu] = isdf(psicr, conj(psicr), ccrank_mu);
cczeta_mu = F*ccrzeta_mu;
psicc_mu = prod_states(psicr(ccind_mu, :), conj(psicr(ccind_mu, :)));

VA = psivc_mu'*vczeta_mu'*Dcoul*vczeta_mu*psivc_mu/(4.0*pi);
VB = psivc_mu'*vczeta_mu'*Dcoul*conj(vczeta_mu*psivc_mu)/(4.0*pi);

WA0 = psicc_mu'*cczeta_mu'*W*vvzeta_mu*psivv_mu/(4.0*pi);
WB0 = psivc_mu'*vczeta_mu'*W*vczeta_mu*psivc_mu/(4.0*pi);

%WA0head = psicc(1,:)'*W(1,1)*psivv(1,:)/(4.0*pi);
%WB0head = psivc(1,:)'*W(1,1)*psivc(1,:)/(4.0*pi);

%WA0wing1 = psicc(1,:)'*W(1,2:end)*psivv(2:end,:)/(4.0*pi);
%WB0wing1 = psivc(1,:)'*W(1,2:end)*psivc(2:end,:)/(4.0*pi);

%WA0wing2 = psicc(2:end,:)'*W(2:end,1)*psivv(1,:)/(4.0*pi);
%WB0wing2 = psivc(2:end,:)'*W(2:end,1)*psivc(1,:)/(4.0*pi);

%WA0body = psicc(2:end,:)'*W(2:end,2:end)*psivv(2:end,:)/(4.0*pi);
%WB0body = psivc(2:end,:)'*W(2:end,2:end)*psivc(2:end,:)/(4.0*pi);


for iv = 1:nvbands,
  for jc = 1:ncbands,
    for ivp = 1:nvbands,
      for jcp = 1:ncbands,
        rl = (iv-1)*ncbands+jc;
        cl = (ivp-1)*ncbands+jcp;
        rr = (jc-1)*ncbands+jcp;
        cr = (iv-1)*nvbands+ivp;
        WA(rl, cl) = WA0(rr, cr);
%        WAhead(rl, cl) = WA0head(rr, cr);
%        WAwing1(rl, cl) = WA0wing1(rr, cr);
%        WAwing2(rl, cl) = WA0wing2(rr, cr);
%        WAbody(rl, cl) = WA0body(rr, cr);
        rl = (iv-1)*ncbands+jc;
        cl = (ivp-1)*ncbands+jcp;
        rr = (ivp-1)*ncbands+jc;
        cr = (iv-1)*ncbands+jcp;
        WB(rl, cl) = WB0(rr, cr);
%        WBhead(rl, cl) = WB0head(rr, cr);
%        WBwing1(rl, cl) = WB0wing1(rr, cr);
%        WBwing2(rl, cl) = WB0wing2(rr, cr);
%        WBbody(rl, cl) = WB0body(rr, cr);
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

