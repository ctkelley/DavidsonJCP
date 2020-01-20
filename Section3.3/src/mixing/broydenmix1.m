function [vnew,df,dv,cdf] = broydenmix1(vin, vout, beta, df, dv, cdf, iter, ...
                                        mixdim, brank)
%
% usage: [vnew,df,dv,cdf] = broydenmix1(vin, vout, beta, df, dv, iter, ...
%                                       mixdim, brank,cdf);
%
% purpose: Broyden mixing for accelerating the SCF iteration
%
% arguements:
%  input: 
%     vin (3-D array) --- input potential
%    vout (3-D array) --- output potential produced from eigenvalue 
%                         calculation
%    df   (vector)    --- work array that stores the changes in function
%                         evaluations.  The dimension of df is n1*n2*n3 by mixdim,
%                         where [n1,n2,n3]=size(vin).
%    dv   (matrix)    --- work array that stores the changes in input 
%                         potentials.  Its dimension is the same as that of df.
%    cdf  (matrix)    --- work array that stores the changes in input 
%                         potentials.  Its dimension is n123*brank by mixdim.
%    iter (integer)   --- current SCF iteration number
%    mixdim(integer)  --- the number of columns in df and dv. This is the 
%                         the mixing memory. Only information from the
%                         previous mixdim iterations are used.
%    brank (integer)  --- the rank of the Broyden update
%  output:
%    vnew (3-D array) --- new potential to be used for the next SCF iteration
%    df   (vector)    --- updated work array that keeps the changes in 
%                         function evaluations
%    dv   (matrix)    --- updated work array that keeps the changes in 
%                           input potentials
%    cdf  (matrix)    --- work array that stores the changes in input 
%                         potentials.  Its dimension is n123*brank by mixdim.
%  
if (brank > mixdim)
    fprintf('brank = %d, mixdim = %d\n', brank, mixdim);
error('The rank of the Broyden update must be larger than the mixing memory size');
end
[n1,n2,n3]=size(vin);
n123 = n1*n2*n3;
%
% Note that the sign convention used in this function is different
% from that used in other mixing functions!!!
%
%
% function evaluation overwrites vout
%
vout = vin - vout;
bragg = zeros(n123,1);
vupd  = beta*reshape(vout,n123,1);
%
if (iter > 1)
   % compute the changes in function evaluations and steps (changes 
   % in potentials)
   if (iter < mixdim+1)
      dv(:,iter-1) = -dv(:,iter-1) + reshape(vin, n123,1);
   else
      dv(:,mixdim) = -dv(:,mixdim) + reshape(vin,n123,1);
   end
   %
   % the difference between the current and the previous function value
   %
   df = -df + reshape(vout,n123,1);
   %
   % Construct C*df
   %
   yvec = beta*df;
   if (iter > brank+1)
      fprintf('Rank-%d update...\n',brank);
      for jbr = 1:iter-brank-1
         jsel = jbr:jbr+brank-1;
         dvmcdf  = dv(:,jsel) - reshape(cdf(:,jbr),n123,brank);
         G = dv(:,jsel)'*reshape(cdf(:,jbr),n123,brank);
         yvec = yvec + dvmcdf*(G\(dv(:,jsel)'*yvec));
      end
      %
      if (brank > 1)
         cdf0 = reshape(cdf(:,jbr),n123,brank);
         cdf0 = cdf0(:,2:brank);
         cdf0 = cdf0 + dvmcdf*(G\(dv(:,jsel)'*cdf0));
         cdf(:,jbr+1) = reshape([cdf0 yvec],n123*brank,1);
      else
         cdf(:,jbr+1) = yvec;
      end
      %
      % Compute the potential update
      %
      for jbr = 1:iter-brank
         jsel = jbr:jbr+brank-1;
         dvmcdf  = dv(:,jsel) - reshape(cdf(:,jbr),n123,brank);
         G = dv(:,jsel)'*reshape(cdf(:,jbr),n123,brank);
         vupd = vupd + dvmcdf*(G\(dv(:,jsel)'*vupd));
      end 
   else
      %
      % perform rank-1 update until there are enough potential vectors to mix
      %
      fprintf('Rank-1 update...\n');
      %
      % Save all initial df vectors in the first column of cdf
      %
      for jbr = 1:iter-2
         ibeg = (jbr-1)*n123+1;
         iend = jbr*n123;
         dvmcdf = dv(:,jbr)-cdf(ibeg:iend,1);
         G = dv(:,jbr)'*cdf(ibeg:iend,1);
         yvec = yvec + dvmcdf*(G\(dv(:,jbr)'*yvec));
      end
      ibeg = n123*(iter-2)+1;
      iend = n123*(iter-1);
      cdf(ibeg:iend,1) = yvec;
      for jbr = 1:iter-1
         ibeg = (jbr-1)*n123+1;
         iend = jbr*n123;
         dvmcdf = dv(:,jbr)-cdf(ibeg:iend,1);
         G = dv(:,jbr)'*cdf(ibeg:iend,1);
         vupd = vupd + dvmcdf*(G\(dv(:,jbr)'*vupd));
      end
   end
end

if (iter < mixdim+1)
   dv(:,iter) = reshape(vin,n123,1);
else
   % delete the first column of dv
   dv(:,1)=[]; 
   % append a new column at the end
   dv(:,mixdim) = reshape(vin,n123,1);
end
df = reshape(vout,n123,1);
vnew = vin - reshape(vupd,n1,n2,n3);
