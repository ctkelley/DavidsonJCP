function X = ifft3(Xfft)
% WAVEFUN/IFFT3 IFFT3 function for wave function class
%    X = IFFT3(Xfft) returns the inverse fast Fourier transform of the wave
%    function.
%
%    See also Wavefun.

%  Copyright (c) 2016-2017 Yingzhou Li and Chao Yang,
%                          Stanford University and Lawrence Berkeley
%                          National Laboratory
%  This file is distributed under the terms of the MIT License.

fpsi = zeros(Xfft.n1,Xfft.n2,Xfft.n3,ncols(Xfft));
N = Xfft.n1*Xfft.n2*Xfft.n3;
if Xfft.iscompact
    for it = 1:ncols(Xfft)
        fpsi(Xfft.idxnz + (it-1)*N) = Xfft.psi(:,it);
    end
else
    fpsi(:) = Xfft.psi(:);
end

xpsi = ifft3(fpsi);
X = Wavefun(xpsi);
X.trans = Xfft.trans;

end