function [vnew,dvmat,vmat] = pulaymix(vin, vout, beta, dvmat, vmat, iter, mixdim)
%
% usage: [vnew,dvmat,vmat] = pulaymix(vin, vout, dvmat, vmat, iter, mixdim);
%
% purpose: Pulay mixing for accelerating the SCF iteration
%
% arguements:
%  input: 
%     vin (3-D array) --- input potential
%    vout (3-D array) --- output potential produced from eigenvalue 
%                         calculation
%    dvmat (matrix)   --- work array that stores the previous vout - vin.
%                         Its dimension is the same as that of vmat.
%    vmat (matrix)    --- work array that stores the previous input potentials
%                         The dimension of vmat is n1*n2*n3 by mixdim,
%                         where [n1,n2,n3]=size(vin).
%    iter (integer)   --- current SCF iteration number
%    mixdim(integer)  --- the number of columns in dvmat and vmat. This is the 
%                         the mixing memory. Only information from the
%                         previous mixdim iterations are used.
%  output:
%    vnew (3-D array) --- new potential to be used for the next SCF iteration
%    dvmat   (matrix) --- updated work array that keeps the previous
%                         vout-vin
%    vmat   (matrix)  --- updated work array that keeps the previous 
%                         input potentials
%  
fprintf('\n Using Pulay mixer... \n');
[n1,n2,n3]=size(vin);
n123 = n1*n2*n3;
%
% construct and modify vmat and dvmat
%
if (iter <= mixdim)
   vmat(:,iter) = reshape(vin,n123,1);
   dvmat(:,iter) = reshape(vout - vin,n123,1);
else
   % delete the first leading column
   vmat(:,1) = []; dvmat(:,1) = [];
   % append vin and vout-vin at the end of vmat and dvmat;
   vmat(:,mixdim) = reshape(vin,n123,1);
   dvmat(:,mixdim) = reshape(vout-vin,n123,1);   
end

% thresh = 1e-6;        % seems to cause crashes in MD
thresh = 1.0e-14; 
if (iter > 1)
   ibeg = 2;
   iend = min(iter,mixdim);
   A = dvmat(:,ibeg:iend);
   b = dvmat(:,ibeg-1);
   ncols = size(A,2);
   A = repmat(b,1,ncols) - A;
   [U,S,V]=svd(A,0);
   while ( S(ncols,ncols) <  thresh*S(1,1) )
      ibeg = ibeg + 1;
      if (ibeg > iend) 
         fprintf('warning: no mixing...\n'); 
         vnew = vout;
         return;
      end
      A = dvmat(:,ibeg:iend);
      b = dvmat(:,ibeg-1);
      ncols = size(A,2);
      A = repmat(b,1,ncols) - A;
      [U,S,V]=svd(A,0);
   end
   g2 = V*(S\(U'*b));
   g1 = 1-sum(g2);
   g = [g1;g2];
   vnew = reshape((vmat(:,ibeg-1:iend)+beta*dvmat(:,ibeg-1:iend))*g,n1,n2,n3); 
else
   vnew = vout;
end
