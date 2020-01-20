function [X,lambda] = diagbyeigs(mol,H,nocc,tol,maxitr)
%
%   [X,lambda] = diagbyeigs(mol,H,nocc,tol,maxitr) Sovle the linear 
%   eigenvalue problem HX = X*Lambda approximate by the MATLAB's eigs 
%   function
%
% Input:
%   mol    a Molecule object
%   H      the Hamiltonian object associated with mol
%   nocc   number of occupied states (desired eigenvalues)
%   tol    convergence tolerance
%   maxitr maximum number of restarts allowed in eigs
%
n1 = mol.n1;
n2 = mol.n2;
n3 = mol.n3;
idxnz = H.idxnz;
nnz = length(idxnz);

eigsopts.isreal = false;
eigsopts.maxit  = maxitr;
eigsopts.tol    = tol;
[V,D,flag]=eigs(@(x)eigsmult(mol,H,x),nnz,nocc,'SR',eigsopts);
d = real(diag(D));
[sd,id]=sort(d);
V = V(:,id);
if (flag~=0)
   fprintf('Convergence not reached in eigs!, pause...\');
   pause;
end
X = Wavefun(V,n1,n2,n3,idxnz);
lambda = sd;

%-----------------------------------------------------
function y = eigsmult(mol,H,x)

n1 = mol.n1;
n2 = mol.n2;
n3 = mol.n3;
idxnz = H.idxnz;
X = Wavefun(x,n1,n2,n3,idxnz);
Y = H*X;
y = Y.psi;