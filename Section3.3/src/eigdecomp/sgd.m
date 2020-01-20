function [X,lambda] = sgd(H, S, X, shift, prec, sm, tol, maxit, verbose)
% SGD for solving biggest or smallest eigenpairs problem
%   [X,lambda] = sgd(H, S, X, shift, prec, sm, tol, maxit, verbose)
%   computes the invariant subspace associated with the k smallest eigenvalues 
%   of the Hamiltonian H by solving the optimization problem: min
%   ||-H+xx*|| in F-norm. The optimization can be solved by gradient descent
%   Note that the objective function is a quartic function with respect to
%   X, so the step size can be determinied by exact line search which
%   corresponding to solving a cubic equation.
%
%   See also omm,ommgd,scg
%
%   SGD adopts optimization to solve engenpair problems. The optimization
%   has the form
%   ..math::`\min_{x\in R^{n\times k}}\left \| H-xx^* \right \|`
%   The optimization have a set of global minimizer :math:`x\in span\{v_1,v_2\cdots v_k\}`
%   where v are the eigenvectors corresbonding to the k-biggest eigenvalues
%   if these eigenvalues are positive, meanwhile, this optimization problem
%   doesn't have any other local minima which guarantee the convergence
%   when using gradient method. For smallest eigenvalue problem, just set H
%   to be -H if the smallest eigenvalues are negetive. The step size in
%   optimization can also be determined by exact line search, the same as
%   OMM. After obtaining vectors x in the invarinat subspace, by solving
%   eigenpair problem for :math:`x^*x` using MATLAB eig function we are
%   able to solve eigenpair priblem for H. This step only requires small
%   amount of calculation cause we convert a huge matrix to a much smaller
%   one.
%
% Input:
%    H     --- a Hamiltonian object
%    X     --- the initial guess of the wave function in Fourier space
%              (Wavefun object)
%    prec  --- preconditioner
%    sm    --- choose smallest or largest eigenvalue
%    tol   --- tolarence on the relative residual.
%    maxit --- maximum number of iterations allowed.
%    verbose - display info
%
% Output
%    X      --- eigenvectors corresponding to the smallest eigenvalues of H
%               (Wavefun object)
%    lambda --- approximations to the smallest eigenvalues of H

if strcmpi(sm,'largestreal') || strcmpi(sm,'lr') || isempty(sm)
    sgn = 1;
elseif strcmpi(sm,'smallestreal') || strcmpi(sm,'sr')
    sgn = -1;
else
    error('SGD: Input variable sm not accepted.');
end

if ~isempty(S)
    SX = S*X;
else
    SX = X;
end
HX = H*X - shift*SX;

for it = 1:maxit
    XSX = X'*SX;
    G = -sgn*4*HX + 4*SX*XSX;
    for j = 1:size(G,2)
        G(:,j) = prec.*G(:,j);
    end
    if ~isempty(S)
        SG = S*G;
    else
        SG = G;
    end
    HG = H*G - shift*SG;
    XSG = X'*SG;
    GHG = G'*HG;
    GSG = G'*SG;
    
    % Solve ax^3+bx^2+cx+d=0
    d = real(trace((XSX)*(XSG) - sgn*G'*HX));
    c = real(trace((XSG)*(XSG) + XSG*XSG' + XSX*GSG - sgn*GHG));
    b = real(3*trace(XSG*GSG));
    a = real(trace(GSG*GSG));
    gamma = linsearchSGD(a,b,c,d);
    
    if (verbose == 1)
        fprintf('SGD iter %3d, err = %11.3e\n', ...
            it, abs(trace(GSG)));
    end
    
    if (abs(trace(GSG))<tol)
        break;
    end
    
    X = X + gamma*G;
    HX = HX + gamma*HG;
    SX = SX + gamma*SG;
end

XSX = X'*SX;
XSX = (XSX+XSX')/2;
[Q,lambda] = eig(XSX);
[lambda,idx] = sort(diag(lambda),'descend','ComparisonMethod','abs');
Q = Q(:,idx);
X = X*(Q/diag(sqrt(lambda)));
lambda = sgn*lambda+shift;

    function gamma = linsearchSGD(a,b,c,d)
        p = (3*c/a-(b/a)^2)/3;
        q = (2*(b/a)^3-9*b/a*c/a+27*d/a)/27;
        p3 = p/3.0;
        q2 = q/2.0;
        
        delta = p3*p3*p3 + q2*q2;
        if (delta >= 0.0)
            qrtd = sqrt(delta);
            gamma = nthroot(-q2+qrtd,3) + nthroot(-q2-qrtd,3);
        else
            qrtd = sqrt(-delta);
            if (q2 >= 0.0)
                gamma = 2 * sqrt(-p3) * cos((atan2( -qrtd, -q2)-2*pi) / 3.0);
            else
                gamma = 2 * sqrt(-p3) * cos(atan2( qrtd, -q2) / 3.0);
            end
        end
        gamma = gamma - b/a/3;
    end

end
