function [cvg,resfro] = reportconverge(H,X,iter,maxiter,curerr,tol,verbose)
% REPORTCONVERGE report convergence.
%    [cvg,resfro] = REPORTCONVERGE(H,X,iter,maxiter,curerr,tol,verbose)
%    reports the convergence of iterations as scf, dcm, etc. And the
%    residue of Frobenius norm. cvg is the indicator for convergence.
%
%    [cvg,resfro] = REPORTCONVERGE(BH,BX,iter,maxiter,curerr,tol,verbose)
%    reports the convergence of iterations as scf, dcm, etc. And the
%    residue of Frobenius norm. The input is with k-points.
%
%   See also scf, scfcry, dcm.

%  Copyright (c) 2016-2017 Yingzhou Li and Chao Yang,
%                          Stanford University and Lawrence Berkeley
%                          National Laboratory
%  This file is distributed under the terms of the MIT License.

if nargin == 6
    verbose = 0;
end

cvg = false;
calcflag = false;
resfro = -1;

if curerr < tol
    fprintf('Convergence is reached!\n');
    calcflag = true;
    cvg = true;
elseif iter == maxiter
    calcflag = true;
end

if calcflag
    if isa(H,'BlochHam')
        resfro = 0;
        nkpts = H.nkpts;
        for it = 1:nkpts
            Hsub = H{it};
            Xsub = X{it};
            HX = Hsub*Xsub;
            G = Xsub'*HX;
            G = (G+G')/2;
            RES = HX-Xsub*G;
            RES = RES.psi;
            if (verbose==1)
                for j = 1:size(RES,2)
                    fprintf('resnrm = %11.3e\n', norm(RES(:,j)));
                end
                fprintf('---------------------\n');
            end
            resfro = resfro + norm(RES,'fro');
        end
    else
        HX = H*X;
        G = X'*HX;
        G = (G+G')/2;
        RES = HX-X*G;
        RES = RES.psi;
        if (verbose==1)
            for j = 1:size(RES,2)
                fprintf('resnrm = %11.3e\n', norm(RES(:,j)));
            end
            fprintf('---------------------\n');
        end
    end
    resfro = norm(RES,'fro');
end

end