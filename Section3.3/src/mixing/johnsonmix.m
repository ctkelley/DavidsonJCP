function [vnew,df,dv] = johnsonmix(vin, vout, df, dv, iter, mixdim, brank)
%
% usage: [vnew,df,dv] = johnsonmix(vin, vout, df, dv, iter, mixdim, brank);
%
% purpose: Johnson mixing for accelerating the SCF iteration
%
% arguements:
%  input: 
%     vin (3-D array) --- input potential
%    vout (3-D array) --- output potential produced from eigenvalue 
%                         calculation
%    df   (matrix)    --- work array that stores the changes in function
%                         evaluations.  The dimension of df is n1*n2*n3 by mixdim,
%                         where [n1,n2,n3]=size(vin).
%    dv   (matrix)    --- work array that stores the changes in input 
%                         potentials.  Its dimension is the same as that of df.
%    iter (integer)   --- current SCF iteration number
%    mixdim(integer)  --- the number of columns in df and dv. This is the 
%                         the mixing memory. Only information from the
%                         previous mixdim iterations are used.
%    brank (integer)  --- the rank of the Johnson update
%  output:
%    vnew (3-D array) --- new potential to be used for the next SCF iteration
%    df   (matrix)    --- updated work array that keeps the changes in 
%                         function evaluations
%    dv   (matrix)    --- updated work array that keeps the changes in 
%                           input potentials
%  
if (brank > mixdim)
    fprintf('brank = %d, mixdim = %d\n', brank, mixdim);
error('The rank of the Johnson update must be larger than the mixing memory size');
end
[n1,n2,n3]=size(vin);
n123 = n1*n2*n3;
w0 = 0.02;
w1 = 1.0;
%
% Note that the sign convention used in this function is different
% from that used in other mixing functions!!!
%
%
% function evaluation overwrites vout
%
vout = vin - vout;
bragg = zeros(n123,1);
fnew  = reshape(vout,n123,1);
%
if (iter > 1)
   % compute the changes in function evaluations and steps (changes 
   % in potentials)
   if (iter < mixdim+1)
       df(:,iter-1) = (-df(:,iter-1) + reshape(vout,n123,1))*w1/w0;
       dv(:,iter-1) = (-dv(:,iter-1) + reshape(vin, n123,1))*w1/w0;
   else
       df(:,mixdim) = (-df(:,mixdim) + reshape(vout,n123,1))*w1/w0;
       dv(:,mixdim) = (-dv(:,mixdim) + reshape(vin, n123,1))*w1/w0;
   end
end
%
if (iter > 1)
   if (iter < mixdim+1)
      jlast = iter-1;
   else
      jlast = mixdim;
   end 
   %
   for jbr = jlast:-1:brank
      %
      % select the columns used in the update
      %
      colsel = jbr-brank+1:jbr;
      %
      % bragg accumulates the Johnson approximation to Jocobian*vout
      %
      cmat = eye(length(colsel))+df(:,colsel)'*df(:,colsel);
      gamma = cmat\(df(:,colsel)'*fnew);
      bragg = bragg - dv(:,colsel)*gamma;
      %
      fnew = fnew - df(:,colsel)*gamma;
   end
end

if (iter < mixdim+1)
   df(:,iter) = reshape(vout,n123,1);
   dv(:,iter) = reshape(vin,n123,1);
else
   % delete the first column from df and dv
   df(:,1)=[];
   dv(:,1)=[]; 
   % append a new column at the end
   df(:,mixdim) = reshape(vout,n123,1);
   dv(:,mixdim) = reshape(vin,n123,1);
end
vnew = vin - reshape(fnew-bragg,n1,n2,n3);
