function Y = power(X,p)
% WAVEFUN/POWER Power function for wave function class
%    Y = POWER(X,p) returns the wave function of X.^p.
%
%    See also Wavefun.

%  Copyright (c) 2016-2017 Yingzhou Li and Chao Yang,
%                          Stanford University and Lawrence Berkeley
%                          National Laboratory
%  This file is distributed under the terms of the MIT License.

Y = X;
Y.psi = X.psi.^p;

end