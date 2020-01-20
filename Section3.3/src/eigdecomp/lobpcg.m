function [X,lambda,LVEC,RVEC] = lobpcg(H, X0, prec, tol, maxit, verbose)
%
% Usage: [X,lambda,LVEC,RVEC] = lobpcg(H, X0, prec, tol, maxit, verbose);
%
% Purpose:
%    Compute the smallest eigenvalue and eigenvector of the Kohn-Sham
%    Hamiltonian using Knyazav's locally optimal preconditioned conjugate
%    gradient (LOBPCG) algorithm.
%
%    Incremental deflation against converged eigenvectors
%
% Input:
%    H     --- a Hamiltonian object
%    X0    --- the initial guess of the wave function in Fourier space
%              (Wavefun object)
%    prec  --- preconditioner
%    tol   --- tolarence on the relative residual.
%    maxit --- maximum number of iterations allowed.
%
% Output:
%    X      --- eigenvectors corresponding to the smallest eigenvalues of H
%               (Wavefun object)
%    lambda --- approximations to the smallest eigenvalues of H
%    LVEC   --- eigenvalue history (matrix of size m by ncols, where m is the 
%               total number of iterations and ncols is the number of eigenvalues
%               computed.)
%    RVEC   --- residual norm history (matrix of size m by ncols)
%

X      = [];
lambda = [];
LVEC   = [];
RVEC   = [];
%
% get size info
%
ncol = ncols(X0); % number of wavefunctions;
if (ncol <=0) 
   fprintf('lobpcg requires at least one initial wave function!\n');
   return;
end
%
% orthonormalize the initial wave functions.
%
[X,~]=qr(X0,0); 
clear X0;
HX = H*X;
%
P  = []; 
HP = []; 
%
nconv = 0;
%
iter = 1;
resnrm = ones(ncol,1);  % initialize residual norm
%
% --- M A I N    L O O P ---
%  
while (iter <= maxit && nconv < ncol)
  % Rayleigh quotient (approximate eigenvalue, obj func)
  S = X'*HX;
  lambda = eig(S);
  lambda = sort(real(lambda)); 
  R = HX - X*S;
  if (verbose == 1)
     fprintf('LOBPCG iter = %3d\n', iter);
  end 
  %
  % Check for convergence
  %
  for j = 1:ncol 
     resnrm(j) = norm(R(:,j));
     if (verbose == 1)
        fprintf('eigval(%2d) = %11.3e, resnrm = %11.3e\n', j, lambda(j), resnrm(j));
     end
  end
  iconv = find(abs(resnrm)<=tol);
  %
  % lock Ritz vectors that satisfy a more stringent convergence tolerance
  %
  tfudge = 1e10;
  ilock = find(abs(resnrm) <= tol/tfudge);
  iact  = find(abs(resnrm) > tol/tfudge);
  %
  % nconv is generally larger than nlock
  % 
  nconv = length(iconv);
  nlock = length(ilock);
  nact  = length(iact);
  %fprintf('nlock = %d, nconv = %d\n', nlock,nconv);
  %
  % test for convergence
  %
  if (nconv >= ncol) 
     break; 
  end 
  % 
  LVEC(iter,1:ncol) = lambda.';
  RVEC(iter,1:ncol) = abs(resnrm)';
  %
  % apply the preconditioner prec
  %
  W = R; 
  for j = 1:ncol
     W(:,j) = prec.*R(:,j);
  end
  if (nlock>0)
     %
     % apply delated preconditioner
     %
     Ylock=X(:,ilock); 
     for j = 1:nlock
        Ylock(:,j) = prec.*Ylock(:,j);
     end
     Xlock=X(:,ilock); 
     Slock = Xlock'*Ylock;
     for j = 1:nact
        W(:,iact(j)) = W(:,iact(j)) - Ylock*(Slock\(Xlock'*W(:,iact(j))));
     end
     %
     % project the deflated Hamiltonian (I-Xlock*Xlock')*H*(I-Xlock*Xlock')
     %
     X(:,iact) = X(:,iact) - Xlock*(Xlock'*X(:,iact));
     W(:,iact) = W(:,iact) - Xlock*(Xlock'*W(:,iact));
     HX(:,iact) = H*X(:,iact);
     HW = H*W(:,iact);
     HX(:,iact) = HX(:,iact) - Xlock*(Xlock'*HX(:,iact));
     HW = HW - Xlock*(Xlock'*HW);
     %
     C = W(:,iact)'*W(:,iact); C = (C+C')/2;
     R = chol(C);
     W(:,iact) = W(:,iact)/R;  
     HW = HW/R; 
     %
     Q  = [X(:,iact) W(:,iact)];
     HQ = [HX(:,iact) HW];
  else
     HW = H*W;
     %
     C  = W'*W; C = (C+C')/2;
     R  = chol(C);
     W  = W/R;  
     HW = HW/R; 
     %
     Q  = [X W];
     HQ = [HX HW];
     if (iter > 1)
        Q  = [Q P];
        HQ = [HQ HP];
     end 
  end
  T = Q'*(HQ); T = (T+T')/2;
  G = Q'*Q; G = (G+G')/2;
  condG = cond(G);
%  if (verbose) 
%     fprintf('cond(G) = %11.3e\n', condG);
%  end
  if (condG > 1.0e12) 
     fprintf('cond(G) = %11.3e\n', condG);
     break;
  end
  [S,D] = eig(T, G,'chol');
  [~,id]= sort(real(diag(D)));
  if (nlock > 0)
     U = S(:,id(1:nact));
     Xact = Q*U;
     X = [Xlock Xact];
     T = X'*(H*X); T = (T+T')/2;
     [S,D] = eig(T);
     [~,id]= sort(real(diag(D)));
     X = X*S(:,id);
     HX = H*X;
  else
     U = S(:,id(1:ncol));
     X = Q*U;
     HX = HQ*U; 
     if (iter > 1)
       set2 = ncol+1:2*ncol;
       set3 = 2*ncol+1:3*ncol; 
       P  = W*U(set2,:)  + P*U(set3,:);
       HP = HW*U(set2,:) + HP*U(set3,:);
       %
%       C = P(:,iact)'*P(:,iact); C = (C + C')/2;
       C = P'*P; C = (C + C')/2;
       R = chol(C);
       P = P/R;
       HP = HP/R;
     else
       P  = W;
       HP = HW;
     end
  end 
  %
  iter = iter + 1;
  if (verbose == 1)
     fprintf('\n'); 
  end
  %pause
end
S = X'*HX; S = (S+S')/2;
[Q,D] = eig(S);
[lambda,id] = sort(real(diag(D)));
X = X*Q(:,id);
