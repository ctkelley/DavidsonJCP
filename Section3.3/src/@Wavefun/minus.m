function Z = minus(X,Y)
% WAVEFUN/MINUS Minus function for wave function class
%    Z = MINUS(X,Y) returns the wave function of X - Y.
%
%    See also Wavefun.

%  Copyright (c) 2016-2017 Yingzhou Li and Chao Yang,
%                          Stanford University and Lawrence Berkeley
%                          National Laboratory
%  This file is distributed under the terms of the MIT License.

if isa(X,'Wavefun') && isa(Y,'Wavefun')
    Z = X;
    Z.psi = X.psi - Y.psi;
    return;
end

if isa(X,'Wavefun')
    Z = X;
    Z.psi = X.psi - Y;
    return;
end

if isa(Y,'Wavefun')
    Z = Y;
    Z.psi = X - Y.psi;
    return;
end

end