function [X1, X2, Lambda] = BSE_complex(A, B)
% The function BSE_complex() computes the structured spectral decomposition
% of a complex Bethe--Salpeter Hamiltonian:
%
%     [A, B; -conj(B), -conj(A)]
%     = [X1, conj(X2); X2, conj(X1)]
%       * blkdiag(Lambda, -Lambda)
%       * [X1, -conj(X2); -X2, conj(X1)]',
%
% where Lambda is diagonal and positive definite.
% The matrix [A, B; conj(B), conj(A)] is required to be positive definite.

i = sqrt(-1);
n = size(A, 1);
I = eye(n);
Q = [I, -i*I; I, i*I]/sqrt(2);
At = [real(A), imag(A); -imag(A), real(A)];
Bt = [real(B), -imag(B); -imag(B), -real(B)];
M = At+Bt;
L = chol(M, 'lower');
J = [0*I, I; -I, 0*I];
W = L'*J*L;
W = (W-W')/2;
[U, T] = hess(W);
t = -diag(T, -1);
T = diag(t, 1) + diag(t, -1);
[V, D] = eig(T);
Lambda = D(n+1:end, n+1:end);
X = Q*(L'\U)*diag(i.^(0:2*n-1))*V;
%[V, D] = eig(-i*W);
%Lambda = D(n+1:end, n+1:end);
%X = Q*(L'\V);
X1 = X(1:n, n+1:end)*sqrt(Lambda);
X2 = X(n+1:end, n+1:end)*sqrt(Lambda);

return
