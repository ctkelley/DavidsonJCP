function VX = applyNP(H,X)
%
% Perform the multiplication of the density-dependent part (Hartree + xc)
% of a Hamiltonian (H) and wavefunction(s) X.
%
% usage: VX = applyNP(H,X);
%

idxnz = H.idxnz;
n1 = X.n1;
n2 = X.n2;
n3 = X.n3;
ncol = ncols(X);
Xr = ifft3(X);
VXr = repmat(H.vnp(:),1,ncol).*Xr;
VX = fft3(VXr);
if X.iscompact
    VX = compress(VX,idxnz);
else
    gm = ones(n1*n2*n3,1);
    gm(idxnz) = 0;
    VX(gm,:) = 0;
end

end
