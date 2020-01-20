function [X1, X2, Lambda] = BSE_real(A, B)
% The function BSE_real() computes the structured spectral decomposition
% of a real Bethe--Salpeter Hamiltonian:
%
%     [A, B; -B, -A]
%     = [X1, X2; X2, X1]
%       * blkdiag(Lambda, -Lambda)
%       * [X1, -X2; -X2, X1]',
%
% where Lambda is diagonal and positive definite.
% The matrix [A, B; B, A] is required to be positive definite.

n = size(A, 1);
M = A+B;
K = A-B;
L1 = chol(M, 'lower');
W = L1'*K*L1;
W = (W+W')/2;
[V, D] = eig(W);
Lambda = sqrt(D);
Phi = L1*V;
Psi = L1'\V;
X1 = (Psi*Lambda+Phi)/2/sqrt(Lambda);
X2 = (Psi*Lambda-Phi)/2/sqrt(Lambda);

return
