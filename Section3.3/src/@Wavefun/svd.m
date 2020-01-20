function [U,S,V]=svd(X,varargin)
% WAVEFUN/SVD Singular value decomposition.
%    [U,S,V] = SVD(X) produces a diagonal matrix S, of the same 
%    dimension as X and with nonnegative diagonal elements in
%    decreasing order, and unitary wave function U and unitary matrix V
%    so that X = U*S*V'.
%  
%    S = svd(X) returns a vector containing the singular values.
%  
%    [U,S,V] = svd(X,0) produces the "economy size"
%    decomposition. If X is m-by-n with m > n, then only the
%    first n columns of U are computed and S is n-by-n.
%    For m <= n, svd(X,0) is equivalent to svd(X).
%  
%    [U,S,V] = svd(X,'econ') also produces the "economy size"
%    decomposition. If X is m-by-n with m >= n, then it is
%    equivalent to svd(X,0). For m < n, only the first m columns 
%    of V are computed and S is m-by-m.
%
%    See also Wavefun.

%  Copyright (c) 2016-2017 Yingzhou Li and Chao Yang,
%                          Stanford University and Lawrence Berkeley
%                          National Laboratory
%  This file is distributed under the terms of the MIT License.

if X.trans
    Xpsi = X.psi';
else
    Xpsi = X.psi;
end

[Umat,S,Vmat] = svd(Xpsi,varargin{:});

if X.trans
    U = Umat;
    V = X;
    V.psi = Vmat;
else
    U = X;
    U.psi = Umat;
    V = Vmat;
end

end