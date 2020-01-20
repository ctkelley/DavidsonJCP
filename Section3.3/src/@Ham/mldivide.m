function X = mldivide(H,B)
% mldivide  Overload the backslash operator for Ham class
%    X = H/B returns a wavefun corresponding to H*X = B.
%
%    See also Ham, Wavefun.

%  Copyright (c) 2016-2017 Yingzhou Li and Chao Yang,
%                          Stanford University and Lawrence Berkeley
%                          National Laboratory
%  This file is distributed under the terms of the MIT License.


n1 = H.n1;
n2 = H.n2;
n3 = H.n3;
idxnz = H.idxnz;

% building the different arrays in this case
xArray = zeros(length(idxnz), ncols(B));


% looping over all the right-hand sides
% TODO: this can be done in parallel
for ii=1:ncols(B)
    % solving the problems using GMRES
    xArray(:,ii) = gmres( @(p) mult(H,p), B.psi(:,ii), [], 1e-6, 300 );
    
end

% Building the resulting Wavefield
X = Wavefun(xArray,n1,n2,n3,idxnz);

end

%-Auxiliary function to leverage the matvec multiplication of the Ham class 
function y = mult(H,x)

n1 = H.n1;
n2 = H.n2;
n3 = H.n3;
idxnz = H.idxnz;

% storing x as a 
X = Wavefun(x,n1,n2,n3,idxnz);
Y = H*X;
y = Y.psi;

end