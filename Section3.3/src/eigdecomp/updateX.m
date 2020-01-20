function [X,ev] = updateX(mol, H, X, prec, options)
%
% usage: [X,ev] = updateX(mol, H, X, prec, options);
%
% purupse: update the wavefunctions by computing the invariant
%          subspace associated with the lowest eigenvalues of H
%
eigmethod = options.eigmethod;
verbose   = options.verbose;
if (any(strcmp(options.verbose,{'off';'OFF';'Off'})))
    verbose = 0;
else
    verbose = 1;
end
ncol = ncols(X);

switch lower(eigmethod)
    case {'lobpcg'}
        cgtol = options.cgtol;
        maxcgiter = options.maxcgiter;
        [X, ev, ~, ~] = lobpcg(H, X, prec, cgtol, maxcgiter,verbose);
    case {'eigs'}
        eigstol = options.eigstol;
        maxeigsiter = options.maxeigsiter;
        [X, ev] = diagbyeigs(mol, H, ncol, eigstol, maxeigsiter);
    case {'chebyfilt'}
        v0 = genX0(mol,1);
        degree = options.degree;
        [T,~,~]=lanczos(H,v0,2*ncol);
        d = sort(real(eig(T)));
        lb = d(ncol+1);
        ub = 1.01*d(2*ncol);
        fprintf('lb = %11.3e, ub = %11.3e, degree = %d\n', ...
            lb, ub, degree);
        Y = chebyfilt(H,X,degree,lb,ub);
        [X,~]=qr(Y,0);
    case {'omm'}
        ommtol = options.cgtol;
        maxommiter = options.maxeigsiter;
        % run a few Lanczos iteration to get an upper bound of the spectrum
        v0 = genX0(mol,1);
        maxlan = 10;
        [T,V,~]=lanczos(H,v0,maxlan);
        T = (T+T')/2;
        [S,D]=eig(T);
        z = V*S(:,maxlan);
        eigub = D(maxlan,maxlan);
        eigub = eigub + norm(H*z-eigub*z);
        [X, ev] = omm(H, [], X, eigub, prec, ommtol, maxommiter,verbose);
    case {'ommgd'}
        ommtol = options.cgtol;
        maxommiter = options.maxeigsiter;
        v0 = genX0(mol,1);
        maxlan = 10;
        [T,V,~]=lanczos(H,v0,maxlan);
        T = (T+T')/2;
        [S,D]=eig(T);
        z = V*S(:,maxlan);
        eigub = D(maxlan,maxlan);
        eigub = eigub + norm(H*z-eigub*z);
        [X, ev] = ommgd(H, [], X, eigub, prec, ommtol, maxommiter,verbose);
    case {'sgd'}
        sgdtol = options.cgtol;
        maxsgditer = options.maxeigsiter;
        v0 = genX0(mol,1);
        maxlan = 10;
        [T,~,~]=lanczos(H,v0,maxlan);
        T = (T+T')/2;
        [~,D]=eig(T);
        shift = D(maxlan,maxlan);
        X = X*diag(sqrt(abs(diag(X'*(H*X))-shift)));
        [X, ev] = sgd(H, [], X, shift, prec, 'sr', sgdtol, maxsgditer,verbose);
    case {'scg'}
        scgtol = options.cgtol;
        maxscgiter = options.maxeigsiter;
        v0 = genX0(mol,1);
        maxlan = 10;
        [T,~,~]=lanczos(H,v0,maxlan);
        T = (T+T')/2;
        [~,D]=eig(T);
        shift = D(maxlan,maxlan);
        X = X*diag(sqrt(abs(diag(X'*(H*X))-shift)));
        [X, ev] = scg(H, [], X, shift, prec, 'sr', scgtol, maxscgiter,verbose);
    otherwise
        disp('Unknown method for diagonalizing H! Use eigs');
        [X, ev] = diagbyeigs(mol, H, ncol, eigstol, maxeigsiter);
end
%
if ( verbose ==1 )
    HX = H*X;
    G = X'*HX;
    R = HX-X*G;
    ev = sort(real(eig(G)));
    fprintf('\n');
    ncol = size(X,2);
    resnrm = zeros(1,ncol);
    for j = 1:ncol
        resnrm(j) = norm(R(:,j));
        fprintf('eigval(%2d) = %11.3e, resnrm = %11.3e\n', ...
            j, ev(j), resnrm(j));
    end
end
