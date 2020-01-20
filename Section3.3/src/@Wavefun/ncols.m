function w = ncols(X)
% WAVEFUN/NCOLS Number of columns in the wave function class
%    w = NCOLS(X) returns the number of columns in the wave function.
%
%    See also Wavefun, nrows.

%  Copyright (c) 2016-2017 Yingzhou Li and Chao Yang,
%                          Stanford University and Lawrence Berkeley
%                          National Laboratory
%  This file is distributed under the terms of the MIT License.

w = size(X.psi,2);

end
