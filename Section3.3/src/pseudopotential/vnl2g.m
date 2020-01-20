function [tab,D] = vnl2g(pp,ecut)
% VNL2G converts non-local pseudo potential at irregular radial grid in
%    real-space to uniform radial grid in G-space.
%    [tab,D] = VNL2G(pp,ecut) calculate the Fourier transform of the
%    non-local pseudo potential to G-space at uniform radial grid in the
%    G-space. The input of this function is the pseudo potential and energy
%    cut. And the function returns the Vnl at G-space in tab and the
%    corresponding middle matrix D. The detailed implementation of the
%    Fourier transform are detailed in the documents.
%
%    See also VLOC2G, GETWQ.

%  Copyright (c) 2015-2016 Yingzhou Li and Chao Yang,
%                          Stanford University and Lawrence Berkeley
%                          National Laboratory
%  This file is distributed under the terms of the MIT License.

dq = 0.01; % gap for interpolation
nqxq = round(sqrt(ecut*2*meDef())/dq+4);

nb = pp.nonloc.nbeta;
beta = pp.nonloc.beta;
hbeta = min(max(pp.nonloc.cutoff_radius_index),length(pp.r));
r = pp.r(1:hbeta);
rab = pp.rab(1:hbeta);
sD = pp.nonloc.dij;
lll = pp.nonloc.lll;

lenl = length(lll);
totall = sum(2*lll+1);

tab = zeros(nqxq,nb);

if nb ==0
    D = [];
    return;
end

for it = 1:nb
    besr = bessel(pp.nonloc.lll(it),r*dq*((1:nqxq)-1));
    aux = repmat(beta(1:hbeta,it),1,nqxq).*besr.*repmat(r,1,nqxq);
    tab(:,it) = simpson(hbeta,aux,rab);
end

% TODO: The following loop could be reduced. Since nb and lll are all
% small, it might be unnecessary to recude it.
D = zeros(totall);
iidx = 0;
for iit = 1:lenl
    ill = pp.nonloc.lll(iit);
    for iitm = -ill:ill
    iidx = iidx+1;
    jidx = 0;
        for jit = 1:lenl
            jll = pp.nonloc.lll(jit);
            for jitm = -jll:jll
                jidx = jidx+1;
                if ill==jll && iitm==jitm
                    D(iidx,jidx) = sD(iit,jit);
                end
            end
        end
    end
end

end
