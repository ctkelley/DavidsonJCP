function Y = fft3(X)
% WAVEFUN/FFT3 FFT3 function for wave function class
%    Xfft = FFT3(X) returns the fast Fourier transform of the wave
%    function.
%
%    See also Wavefun.

%  Copyright (c) 2016-2017 Yingzhou Li and Chao Yang,
%                          Stanford University and Lawrence Berkeley
%                          National Laboratory
%  This file is distributed under the terms of the MIT License.

N = X.n1*X.n2*X.n3;
if X.iscompact
    xpsi = zeros(X.n1,X.n2,X.n3,ncols(X));
    for it = 1:ncols(X)
        xpsi(X.idxnz + (it-1)*N) = X.psi(:,it);
    end
else
    xpsi = reshape(X.psi,X.n1,X.n2,X.n3,ncols(X));
end

fpsi = fft3(xpsi);
Y = Wavefun(fpsi);
Y.trans = X.trans;

end