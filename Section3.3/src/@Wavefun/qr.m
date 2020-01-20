function varargout = qr(X,varargin)
% WAVEFUN/QR Orthogonal-triangular decomposition.
%    [Q,R] = QR(A), where A is a wave function, Q is of size nrows by
%    ncols and R is upper triangular matrix of size ncols by ncols.
%
%    [Q,R] = qr(A,0) produces the "economy size" decomposition.
%
%    [Q,R,E] = qr(A) produces unitary Q, upper triangular R and a
%    permutation matrix E so that A*E = Q*R. The column permutation E is
%    chosen so that ABS(DIAG(R)) is decreasing.
%
%    [Q,R,e] = qr(A,'vector') returns the permutation information as a
%    vector instead of a matrix.  That is, e is a row vector such that
%    A(:,e) = Q*R. Similarly, [Q,R,E] = qr(A,'matrix') returns a permutation
%    matrix E. This is the default behavior.
%
%    [Q,R,E] = qr(A,0) produces an "economy size" decomposition in which E
%    is a permutation vector, so that A(:,E) = Q*R.
%
%    See also Wavefun.

%  Copyright (c) 2016-2018 Yingzhou Li and Chao Yang,
%                          Stanford University and Lawrence Berkeley
%                          National Laboratory
%  This file is distributed under the terms of the MIT License.

if X.trans
    Xpsi = X.psi;
    if nargout == 2
        [Q,Rmat] = qr(Xpsi',varargin{:});
    else
        [Q,Rmat,E] = qr(Xpsi',varargin{:});
    end
    
    R = X;
    R.psi = Rmat';
    R.trans = 1;
    varargout{1} = Q;
    varargout{2} = R;
    if nargout == 3
        varargout{3} = E;
    end
else
    Xpsi = X.psi;
    if nargout == 2
        [Qmat,R] = qr(Xpsi,varargin{:});
    else
        [Qmat,R,E] = qr(Xpsi,varargin{:});
    end
    
    Q = X;
    Q.psi = Qmat;
    varargout{1} = Q;
    varargout{2} = R;
    if nargout == 3
        varargout{3} = E;
    end
end

end