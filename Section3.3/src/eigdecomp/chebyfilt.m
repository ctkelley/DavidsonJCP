function  Y = chebyfilt(H, X, degree, lb, ub)
% CHEBYFILT chebyshef filter for engenvalue problem
%   Y = CHEBYFILT(H, X, degree, lb, ub) Apply p(H) to X to obtain 
%   Y = p(H)*X, where p(lambda) is a shifted and scaled Chebyshev polynomial
%   of degree 'degree'.  The polynomial is bounded between -1 and 1 within 
%   [lb,ub].
%   
%   To solve a eigenvalue problem, the traditional algorithm requires
%   calclulation amount of order :math:`O(n^3)` , which is very expensive
%   when applying on a huge matrix with large n, especially when we only
%   requires some of the matrix's eigenvalues such calculation is
%   unnecessary. So, a few kinds of fast algorithms is provided in floder
%   'eigdecomp' including chebyshef filter, OMM, lobpcg, gmpcg, scg,
%   zolotarev function and use eig function in MATLAB. 
%
%   The chebyshef filter and zolotarev function is used to approximate a
%   function used on matrix H to filt the eigenvalues. It's obvious that a
%   matrix function :math:`f(H)` is equal to :math:`Q*f(\Lambda )Q` (with 
%   some constrains on the smothness of f) where H
%   is a postive difinite matrix and :math:`H=Q*\Lambda Q` is the eigen
%   decomposition of H, Q is a orthogonal matrix and :math:`\Lambda` is a
%   diagonal matrix consisting all H's eigenvalues. So, when we only
%   concern some of the eigenvalues, the function can be chosen as a
%   filter, to be more specific, if we are able to konw the interval which
%   contains the eigenvalues we are interested in, such as several smallest
%   or biggest eigenvalues, the function is chosen as 
%   ..math::f(x)=\left\{\begin{matrix} 1\ \ x\in I\\ 0\ else \end{matrix}\right.
%   The function is not smooth at the edges of interval I but by assuming
%   that there is a gap between the edges and the nearest eigenvalue, the
%   equation still holds. Now the only concern becomes how to approximate
%   such :math:`f(x)`, note that the function is applied to a matrix, and
%   most matrix functions are hard to compute, the common choice is
%   polynomial function. We choose chebyshef polynomials as a basis which
%   are a set of orthogonal polynomials to represent our objective
%   function.
%
%   Once we are able to approximate the objective function, we can use the
%   following algorithm to finish solving generalized engenpair problem::
%
%     Input:Hermitian matrix pensil (A,B) with size :math:`n \times n`,
%     spectrum range (a,b) including the eigenvalue we are interested in,
%     number of eigenpairs :math:`n_\lambda` and the approximate filter
%     function :math:`R_ab`
%     Output: diagonal matrix :math:`\Lambda` with diagonal entries being
%     eigenvalues of (A,B) on (a,b), and V the corresbonding eigenvectors
%     Choose random orthonormal vectors :math:`Q \in F^{N\times(n_lambda+k)}`
%     while not convergent do
%        :math:`Y=R_ab(B^{-1}A)`
%        :math:`\widetilde{A}=Y^*AY,\widetilde{B}=Y^*BY`
%        solve :math:`\widetilde{A}\widetilde{Q}=\widetilde{\Lambda}\widetilde{B}\widetilde{Q}` for :math:`\widetilde{\Lambda}` and :math:`\widetilde{Q}`
%        updage :math:`Q=Y\widetilde{Q}`
%     end
%     :math:`\widetilde{I}=\{i\left | a<\widetilde{\Lambda}(i,i)<b \}`
%     :math:`\lambda=\widetilde{\Lambda}(\widetilde{I},\widetilde{I})`
%     :math:`V=Q(:,\widetilde{I})`
%
%   Note that the algorithm oonly requires to sove a eigenvalue problem for
%   a matrix with size :math:`n_\lambda+k` which is much cheaper in
%   situation such as KSSLOV
%
%
% Input:
%   H (Ham object)      Hamiltonian 
%   X (Wavefun object)  Wavefunction
%   degree (integer)    The degree of the Chebyshev polynomial
%   lb (float)          The lower bound of the interval in which the 
%                       polynomial is bounded between -1 and 1.
%   ub (float)          The upper bound of the interval in which the 
%                       polynomial is bounded between -1 and 1.
% Output:  
%
%   Y  (Wavefun object)  Wavefunction Y = p(H)*X;
%
e = (ub-lb)/2;
c = (ub+lb)/2;
sigma = e/(lb-ub);
sigma1 = sigma;
Y = H*X - X*c;
Y = Y*sigma1/e;
for i = 2:degree
   sigma2 = 1/(2/sigma1 - sigma);
   Y1 = (H*Y-Y*c)*2*sigma2/e-X*(sigma*sigma2);
   X = Y;
   Y = Y1;
   sigma = sigma2;
end
