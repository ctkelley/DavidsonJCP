function [X,lambda] = scg(H, S, X, shift, prec, sm, tol, maxit, verbose)
%
% Usage: [X,lambda,LVEC,RVEC] = scg(H, X0, eigub, prec, tol, maxit, verbose)
%
%

if strcmpi(sm,'largestreal') || strcmpi(sm,'lr') || isempty(sm)
    sgn = 1;
elseif strcmpi(sm,'smallestreal') || strcmpi(sm,'sr')
    sgn = -1;
else
    error('SCG: Input variable sm not accepted.');
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

    if it > 1
        beta = sum(sum(conj(G).*(G-G0)))/norm(G0)^2;
        D = G + G0*beta;
        G0 = G;
    else
        D = G;
        G0 = G;
    end
    
    if ~isempty(S)
        SD = S*D;
    else
        SD = D;
    end
    HD = H*D - shift*SD;
    XSD = X'*SD;
    DHD = D'*HD;
    DSD = D'*SD;
    
    % Solve ax^3+bx^2+cx+d=0
    d = real(trace((XSX)*(XSD) - sgn*D'*HX));
    c = real(trace((XSD)*(XSD) + XSD*XSD' + XSX*DSD - sgn*DHD));
    b = real(3*trace(XSD*DSD));
    a = real(trace(DSD*DSD));
    gamma = linsearchSCG(a,b,c,d);
    
    if (verbose == 1)
        fprintf('SGD iter %3d, err = %11.3e\n', ...
            it, abs(trace(DSD)));
    end
    
    if (abs(trace(DSD))<tol)
        break;
    end
    
    X = X + gamma*D;
    HX = HX + gamma*HD;
    SX = SX + gamma*SD;
end

XSX = X'*SX;
XSX = (XSX+XSX')/2;
[Q,lambda] = eig(XSX);
[lambda,idx] = sort(diag(lambda),'descend','ComparisonMethod','abs');
Q = Q(:,idx);
X = X*(Q/diag(sqrt(lambda)));
lambda = sgn*lambda+shift;

    function gamma = linsearchSCG(a,b,c,d)
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
