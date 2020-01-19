% DAVIDSON
% Shameless hack job to get it working.
%
% function [u,lambda,itc,ihist]=davidson(A,M,u,tol,maxit,lambda0)
% input:
% A = coefficient matrix; M = preconditioner
% u = initial iterate for eigenvector
% maxit = what you think it is
% lambda = initial iterate for eigenvalue
%
% output:
% u = eigenvector; lambda = eigenvalue
% itc = iteration counter
% ithist = relative residual history
%
function [u,lambda,itc,ihist]=davidson(A,M,u,tol,maxit,lambda0)
n=length(u);
if nargin == 5
lambda=u'*A*u/(u'*u);
else
lambda=lambda0;
end
r=A*u - lambda *u;
%nr=norm(r)/sqrt(n);
nr=norm(r);
V=[u];
itc=0;
tol=tol*nr;
nr0=nr;
ihist=nr/nr0;
while nr > tol && itc < maxit
%    ML=M - lambda*speye(n,n);
    t=-M\r;
    p=u+t;
    V=gram(V,p);
    AV=V'*A*V;
    [w,lambda]=eigs(AV,1,'SM');
    u=V*w;
    r=A*u-lambda*u;
%    nr=norm(r)/sqrt(n);
    nr=norm(r);
    ihist=[ihist; nr/nr0];
    itc=itc+1;
end
end

function V=gram(V,u)
w1=u - V*(V'*u);
v1=w1/norm(w1);
w2=v1-V*(V'*v1);
v2=w2/norm(w2);
V=[V,v2];
end
