function Y = ctranspose(X)
% WAVEFUN/CTRANSPOSE Ctranspose function for wave function class
%    Y = CTRANSPOSE(X) returns the conjugate transpose value of the wave
%    function.
%
%    See also Wavefun.

%  Copyright (c) 2016-2017 Yingzhou Li and Chao Yang,
%                          Stanford University and Lawrence Berkeley
%                          National Laboratory
%  This file is distributed under the terms of the MIT License.

Y = X;
Y.trans = X.trans == 0;

end