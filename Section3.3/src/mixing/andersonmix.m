function [vnew,df,dv] = andersonmix(vin, vout, beta, df, dv, iter, mixdim)
%
% usage: [vnew,df,dv] = andersonmix(vin, vout, beta, df, dv, iter, mixdim)
%
% purpose: Anderson mixing for accelerating the SCF iteration
%
% Anderson mixing is a special case for Broyden mixing
% The mixing can be expressed by
%
%  vnew = vin - beta* Jinv * r
%
% where r = vout - vin, and Jinv is an approximation to the inverse
% of the Jacobian.
% 
% For Anderson mixing, Jinv is defined to be
%
% Jinv = I + (dv - df)*pinv(df)
%
% This expression can be derived from  min_Jinv || Jinv - I ||_F
%                                      s.t. dv =  Jinv*df  
%
% arguements:
%  input: 
%     vin (3-D array) --- input potential
%    vout (3-D array) --- output potential produced from eigenvalue 
%                         calculation
%    beta (scalar)    --- mixing parameter between 0 and 1
%    df   (matrix)    --- work array that stores the changes in function
%                         evaluations.  The dimension of df is n1*n2*n3 by mixdim,
%                         where [n1,n2,n3]=size(vin).
%    dv   (matrix)    --- work array that stores the changes in input 
%                         potentials.  Its dimension is the same as that of df.
%    iter (integer)   --- current SCF iteration number
%    mixdim(integer)  --- the number of columns in df and dv. This is the 
%                         the mixing memory. Only information from the
%                         previous mixdim iterations are used.
%  output:
%    vnew (3-D array) --- new potential to be used for the next SCF iteration
%    df   (matrix)    --- updated work array that keeps the changes in 
%                         function evaluations
%    dv   (matrix)    --- updated work array that keeps the changes in 
%                           input potentials
%
[n1,n2,n3]=size(vin);
n123 = n1*n2*n3;
%
% function evaluation overwrites vout
%
vout = vout - vin;
%
iterused = min(iter-1,mixdim);
ipos = iter - 1 - floor((iter-2)/mixdim)*mixdim;
%
if (iter > 1)
   % compute the changes in function evaluations and the step (changes 
   % in potentials)
   df(:,ipos) = df(:,ipos) - reshape(vout,n123,1);
   dv(:,ipos) = dv(:,ipos) - reshape(vin, n123,1);
end
%
vinsave  = vin;
voutsave = vout;
%
if (iter > 1)
   %B = (df(:,1:iterused))'*(df(:,1:iterused));
   %B = (B+B')/2;
   % 
   %work = df(:,1:iterused)'*reshape(vout,n123,1);
   %gammas = B\(work);
   %
   % the use of pseudo-inverse is not the most efficient
   % implementation.
   %
   gammas = pinv(df(:,1:iterused))*reshape(vout,n123,1);  
   %
   % the following loop can be replaced by two simple gemvs
   % and some reshape operations.
   % 
   for i = 1:iterused
      vin  = vin  - gammas(i)*reshape(dv(:,i),n1,n2,n3);
      vout = vout - gammas(i)*reshape(df(:,i),n1,n2,n3);
   end
end

inext = iter - floor((iter - 1) / mixdim) * mixdim;
df(:,inext) = reshape(voutsave,n123,1);
dv(:,inext) = reshape(vinsave,n123,1);

vnew = vin + beta*vout;
