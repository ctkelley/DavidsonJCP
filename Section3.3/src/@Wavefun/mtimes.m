function Z = mtimes(X,Y)
% WAVEFUN/MTIMES Mtimes function for wave function class
%    M = MTIMES(X,Y) returns a matrix as the multiplication of two wave
%    functions.
%
%    Z = MTIMES(M,Y) or Z = MTIMES(X,M) returns a wave function as the
%    multiplication of a matrix or scalar with a wave function.
%
%    Note: the left multiplication of a non-scalar matrix with a wave
%    function returns a regular matrix.
%
%    See also Wavefun.

%  Copyright (c) 2016-2017 Yingzhou Li and Chao Yang,
%                          Stanford University and Lawrence Berkeley
%                          National Laboratory
%  This file is distributed under the terms of the MIT License.

if isa(X,'Wavefun') && isa(Y,'Wavefun')
    if X.trans
        Xpsi = X.psi';
    else
        Xpsi = X.psi;
    end
    if Y.trans
        Ypsi = Y.psi';
    else
        Ypsi = Y.psi;
    end
    Z = Xpsi*Ypsi;
    return;
end

if isa(X,'Wavefun')
    if X.trans
        Z = X.psi'*Y;
    else
        Z = X;
        Z.psi = X.psi*Y;
    end
    return;
end

if isa(Y,'Wavefun')
    if Y.trans
        Z = Y;
        Z.psi = Y.psi*X';
    else
        if isscalar(X)
            Z = Y;
            Z.psi = X * Y.psi;
        else
            Z = X*Y.psi;
        end
    end
    return;
end

end
