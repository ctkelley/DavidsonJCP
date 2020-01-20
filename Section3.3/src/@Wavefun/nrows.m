function h = nrows(X)
% WAVEFUN/NROWS Number of rows in the wave function class
%    h = NROWS(X) returns the number of rows in the wave function.
%
%    See also Wavefun, ncols.

%  Copyright (c) 2016-2017 Yingzhou Li and Chao Yang,
%                          Stanford University and Lawrence Berkeley
%                          National Laboratory
%  This file is distributed under the terms of the MIT License.

h = size(X.psi,1);

end