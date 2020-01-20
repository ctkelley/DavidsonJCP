function Z = mrdivide(X,R)
% WAVEFUN/MRDIVIDE Mrdivide function for wave function class
%    Z = MRDIVIDE(X,R) returns the wave function of X/R.
%
%    See also Wavefun.

%  Copyright (c) 2016-2017 Yingzhou Li and Chao Yang,
%                          Stanford University and Lawrence Berkeley
%                          National Laboratory
%  This file is distributed under the terms of the MIT License.

Z = X;
Z.psi = X.psi/R;

end