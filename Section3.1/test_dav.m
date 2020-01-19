% TEST_DAV
% Generates the iteration statistics for section 3.1
% 
% function ihist=test_dav(n)
% Runs the example in Section 3.1 with h=1/(n+1)
%
function ihist=test_dav(n)
if nargin == 0
  n=8000;
end
h=1/(n+1);
x=h:h:1-h; x=x';
e = ones(n,1);
M = spdiags([-e 2*e -e], -1:1, n, n);
M=M/(h*h);
K=intop(n);
A = M + K;
u=x.*(1-x); tol=1.d-7; 
%u=sqrt(n)*u/norm(u);
u=u/norm(u);
maxit=20;
tic
%[u, lambda,itc, ihist]=davidson(A,M,u,tol,maxit,10.0);
[u, lambda,itc, ihist]=davidson(A,M,u,tol,maxit);
davtime=toc;
tic
[ut,lambdat]=eigs(A,1,'SM');
eigstime=toc;
dt=u-(ut'*u)*ut;
[norm(dt,inf),lambda-lambdat,itc]
ut=ut*sqrt(n)/norm(ut);
u=u*sqrt(n)/norm(u);
debug=0;
if debug==1
[davtime,eigstime]
plot(x,u,'-',x,ut,'--')
[lambdat,lambda]
ihist
end
end

function K=intop(n)
h=1/(n+1);
x=h:h:1-h; x=x';
%K=20*cos(-10*pi*abs(x-x'))*h;
K=-exp(.5*abs(x-x'))*h;
end

