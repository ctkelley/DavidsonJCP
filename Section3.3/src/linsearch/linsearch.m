function [XP, idid, lambda]=linsearch(fc,qp0,ft,X0,U,S,V,alpha,maxarm)
%
% Adapted from C. T. Kelley's polyline code
%
% function [XP, idid, lambda]=linsearch(fc,qp0,ft,X0,U,S,V,alpha,maxarm);
%
% polynomial line search, call after first point is rejected
%
% Input: 
%        q0      = current function value
%        qp0     = current directional derivative
%        ft      = trial function (rejected value)
%        X       = current point
%        [U,S,V] = svd(search direction)
%        alpha   = parameter used in total energy function 
%                  evaluation nlham.
%        maxarm  = maximum number of step length reductions   
%
% Output: XP = successful new point (if it exists)
%       idid = number of calls to f (if line search succeeds) or
%              -1 if line search fails.
%
% Requires: polymod.m
%
% line search parameters that everyone uses
%
alp=1.d-4; blow=.1; bhigh=.5;
%
% Set up the search
%
q0=fc; qc=ft; lamc=1; iarm=0; numf=0;
ds = diag(S);
fgoal=q0+alp*lamc*qp0;
while ft > fgoal
    iarm=iarm+1;
    if iarm==1  % quadratic
       lambda=polymod(q0, qp0, lamc, qc, blow, bhigh);
    else
       lambda=polymod(q0, qp0, lamc, qc, blow, bhigh, lamm, qm);
    end
    qm=qc; lamm=lamc; lamc=lambda;
    cs = cos(ds*lambda);
    sn = sin(ds*lambda); 
    X = X0*V*diag(cs)*V'+U*diag(sn)*V';
    [H,rho,ft]=nlham(X,alpha); numf = numf+1; qc=ft;
    if(iarm > maxarm)
         disp(' line search failure'); idid=-1; XP=X0;
    return; end
    fgoal=fc+alp*lamc*qp0;
end
XP=X; idid=numf;
