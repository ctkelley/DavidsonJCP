function X = compress(X,idxnz)
% WAVEFUN/COMPRESS Compress function for wave function class
%    X = COMPRESS(X,idxnz) returns the compact wave function with the
%    non-zero index idxnz.
%
%    See also Wavefun.

%  Copyright (c) 2016-2017 Yingzhou Li and Chao Yang,
%                          Stanford University and Lawrence Berkeley
%                          National Laboratory
%  This file is distributed under the terms of the MIT License.

X.idxnz = idxnz;
X.psi = X.psi(idxnz,:);
X.iscompact = 1;

end