function beta = sum(X,varargin)
% WAVEFUN/SUM Sum function for wave function class
%    beta = SUM(X) returns the sum of the wave function.
%
%    beta = SUM(X,dim) returns the sum of the wave function along the dim
%    dimension.
%
%    See also Wavefun.

%  Copyright (c) 2016-2017 Yingzhou Li and Chao Yang,
%                          Stanford University and Lawrence Berkeley
%                          National Laboratory
%  This file is distributed under the terms of the MIT License.

beta = sum(X.psi,varargin{:});

end