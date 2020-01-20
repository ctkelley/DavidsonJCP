function [C,lambda] = ommgd(H, S, C, shift, prec, tol, maxit, verbose)
%
% Usage: [X,lambda,LVEC,RVEC] = ommdg(H, X0, eigub, prec, tol, maxit, verbose)
%
% Purpose:
%    Compute the invariant subspace associated with the k smallest
%    eigenvalues of the Hamiltonian H using Orbital Minimization,
%    where k is the number of columns in X0
%
% Input:
%    H     --- a Hamiltonian object
%    X0    --- the initial guess of the wave function in Fourier space
%              (Wavefun object)
%    eigub --- upper bound of the spectrum of H
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

HC = H*C - shift*C;
Hw = C'*HC;
if (~isempty(S))
    SC = S*C;
else
    SC = C;
end
Sw = C'*SC;

for it = 1:maxit
    
    G = - 8*HC + 4*SC*Hw + 4*HC*Sw; %Negative gradient
    
    % apply the preconditioner prec
    for j = 1:size(G,2)
        G(:,j) = prec.*G(:,j);
    end
    
    Hw2 = G'*HC;
    HG = H*G - shift*G;
    Hw3 = G'*HG;
    Sw2 = G'*SC;
    if (~isempty(S))
        SG = S*G;
    else
        SG = G;
    end
    Sw3 = G'*SG;
    
    % Solve ax^3+bx^2+cx+d=0
    d = real(8*trace(Hw2) - 4*trace(Sw2*Hw) - 4*trace(Sw*Hw2));
    c = real(8*trace(Hw3) - 4*trace(Sw3*Hw) - 4*trace(Sw*Hw3) ...
        - 8*trace(Sw2*Hw2) - 8*trace(Sw2'*Hw2));
    b = real(-12*trace(Sw3*Hw2) - 12*trace(Sw2*Hw3));
    a = real(-8*trace(Sw3*Hw3));
    gamma = linsearchOMM(a,b,c,d);
    
    % update
    C = C + gamma*G;
    HC = HC + gamma*HG;
    Hw = Hw + gamma*(Hw2 + Hw2') + gamma^2*Hw3;
    SC = SC + gamma*SG;
    Sw = Sw + gamma*(Sw2 + Sw2') + gamma^2*Sw3;
    
    if (verbose == 1)
        fprintf('OMM iter = %3d, rel err = %11.3e\n', it, abs(trace(Sw3)));
    end
    if (abs(trace(Sw3)) < tol)
        break;
    end
    
end

% Rayleigh-Ritz for now (can be removed later)
Hw = (Hw+Hw')/2;
Sw = (Sw+Sw')/2;
[Q,lambda] = eig(Hw,Sw);
[lambda,id] = sort(diag(lambda));
lambda = lambda + shift;
C = C*Q(:,id);

    function gamma = linsearchOMM(a,b,c,d)
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

