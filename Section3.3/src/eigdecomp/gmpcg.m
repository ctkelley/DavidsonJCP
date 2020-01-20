function [X,d] = gmpcg(X0,alpha,maxit);
%
% usage: [X,d] = gmpcg(X0,alpha,maxit);
%
% Solve the toy nonlinear eigenvalue problem using Edelman's
% Grassmann (preconditioned) conjugate gradient
%
tol = 1e-10;
[n,ncols] = size(X0);
[X,R] = qr(X0,0);
%
[H,rho0,etot0] = nlham(X,alp+ha);
fprintf('Total energy = %11.3e\n', etot0);
HX = H*X;
G0 = HX - X*(X'*HX);
P0 = -G0;
X0 = X;
L  = trid(-1,2,-1,n); 
%
for iter = 1:maxit
   [U,S,V]=svd(P0,0);
   ds = diag(S);
   cs = cos(ds);
   sn = sin(ds); 
   %
   % take a full step first
   %
   X = X0*V*diag(cs)*V' + U*diag(sn)*V';
   %
   [H,rho,etot]=nlham(X,alpha);
   tmin = 1;
   %
   % perform a line search if taking a full step 
   % does not lead to a decrease in etot.
   %
   if (etot > etot0)
     display('need line search...');
     maxarm = 10;
     fc = etot0;
     ft = etot;
     qp0 = trace(V*S*(U'*L*X0)) + alpha*rho0'*(L\diag(X0*V*S*U'));
     [X, idid, tmin]=linsearch(fc,qp0,ft,X0,U,S,V,alpha,maxarm);
     fprintf('tmin = %11.3e\n', tmin);
     if (idid < 0) 
        fprintf('Line Search failed...\n');
        break;
     end
     ds = diag(S)*tmin;
     cs = cos(ds);
     sn = sin(ds);
     X = X0*V*diag(cs)*V'+ U*diag(sn)*V';
     [H,rho,etot] = nlham(X,alpha);
   end
   %fprintf('orth = %11.3e\n', norm(X'*X-eye(ncols),'fro'))
   %
   % compute new gradient
   %
   HX = H*X;
   G  = HX - X*(X'*HX);
   resnrm  = norm(G,'fro');
   err     = norm(rho-rho0);
   orth    = norm(X'*X - eye(ncols), 'fro');
   fprintf('iter = %d\n', iter);
   fprintf('total energy   = %11.3e\n', etot);
   fprintf('norm(rho-rho0) = %11.3e\n', err);
   fprintf('resnrm         = %11.3e\n', resnrm);
   fprintf('orth           = %11.3e\n', orth);
   fprintf('----------\n');
   pause;
   if ( err < tol & resnrm < tol )
     break;
   end
   %
   % parallel transport the previous search direction and gradient
   %
   Pt = (-X0*V*diag(sn) + U*diag(cs))*S*V';
   UG = U'*G0;
   Gt = G0 - ( X0*V*diag(sn)*UG + U*(UG - diag(cs)*UG) );
   %
   % compute new search direction
   %
   gam = sum(sum((G-Gt).*G))/sum(sum(G0.*G0));
   P   = -G + gam*Pt;
   %
   rho0 = rho; P0 = P; G0 = G; etot0 = etot; 
   %[X0,R0] = qr(X,0);
   X0  = X;
end
%
[V,D]=eig(H);
d = diag(D(1:ncols,1:ncols));
[H,rho,etot] = nlham(X,alpha);
HX = H*X;
G = HX -X*(X'*HX);
resnrm = norm(G,'fro');
fprintf('Total energy = %11.3e\n', etot);
fprintf('resid norm   = %11.3e\n', norm(G,'fro'));
