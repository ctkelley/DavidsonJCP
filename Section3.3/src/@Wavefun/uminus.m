function Y = uminus(X)
% WAVEFUN/UMINUS Uminus function for wave function class
%    Y = UMINUS(X) returns a wave function as -X.
%
%    See also Wavefun.

%  Copyright (c) 2016-2017 Yingzhou Li and Chao Yang,
%                          Stanford University and Lawrence Berkeley
%                          National Laboratory
%  This file is distributed under the terms of the MIT License.

Y = X;
Y.psi = -X.psi;

end