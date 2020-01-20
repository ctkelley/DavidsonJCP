function beta = norm(X,varargin)
% WAVEFUN/NORM Norm function for wave function class
%    beta = NORM(X,opt) returns the norm of the wave functions.
%
%    See also Wavefun.

%  Copyright (c) 2016-2017 Yingzhou Li and Chao Yang,
%                          Stanford University and Lawrence Berkeley
%                          National Laboratory
%  This file is distributed under the terms of the MIT License.

beta = norm(X.psi,varargin{:});

end