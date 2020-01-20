function Xconj = conj(X)
% WAVEFUN/CONJ Conj function for wave function class
%    Xconj = CONJ(X) returns the conjugate value of the wave function.
%
%    See also Wavefun.

%  Copyright (c) 2016-2017 Yingzhou Li and Chao Yang,
%                          Stanford University and Lawrence Berkeley
%                          National Laboratory
%  This file is distributed under the terms of the MIT License.

Xconj = X;
Xconj.psi = conj(X.psi);

end