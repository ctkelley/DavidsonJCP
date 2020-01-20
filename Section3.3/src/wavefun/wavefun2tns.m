function T = wavefun2tns(X,j)
% WAVEFUN2TNS  Convert a slice of wave function to tensor.
%   T = WAVEFUN2TNS(X) converts the all slices of wave function to tensor.
%
%   T = WAVEFUN2TNS(X,j) converts the jth slice of wave function to tensor.
%
%
%   See also Wavefun.

%  Copyright (c) 2016-2017 Yingzhou Li and Chao Yang,
%                          Stanford University and Lawrence Berkeley
%                          National Laboratory
%  This file is distributed under the terms of the MIT License.

if nargin < 2
    j = 1:ncols(X);
end

n1 = X.n1;
n2 = X.n2;
n3 = X.n3;
Xpsi = X.psi;
if X.iscompact
    T = zeros(n1*n2*n3,numel(j));
    T(X.idxnz,:) = X.psi(:,j);
else
    T = X.psi(:,j);
end

T = reshape(T,n1,n2,n3,[]);

end
