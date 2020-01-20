function Xabs = abs(X)
% WAVEFUN/ABS Abs function for wave function class
%    Xabs = ABS(X) returns the absolute value of the wave function.
%
%    See also Wavefun.

%  Copyright (c) 2016-2017 Yingzhou Li and Chao Yang,
%                          Stanford University and Lawrence Berkeley
%                          National Laboratory
%  This file is distributed under the terms of the MIT License.

Xabs = X;
Xabs.psi = abs(X.psi);

end
