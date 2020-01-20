function y = mtimes(F,x)
%
% Perform the inverse Fourier transform from the reciprocal 
% space to the real space
%
% usage: y = mtimes(F,x);
%
if (nargin == 2)
   idxnz = F.idxnz;
   n1 = F.n1;
   n2 = F.n2;
   n3 = F.n3;
   n123 = n1*n2*n3;
   vol = F.vol;
   if ( isa(x,'numeric') )
      [nrows, ncols ] = size(x);
      if (F.inverse)  
         ng = length(idxnz);
         if ( nrows ~= ng )
           error('the number of rows in x does not match with the KSFFT object, nrows = %d, ng = %d', nrows, ng);
         end
         a3d = zeros(n1,n2,n3);
         y = zeros(n123,ncols);
         for j = 1:ncols 
            a3d(idxnz) = x(:,j);
            f3d = ifftn(a3d);
            y(:,j) = f3d(:);
            a3d(idxnz) = 0;
         end;
         y = y*n123/vol;
      elseif (F.forward)
         if ( nrows ~= n123 )
           error('the number of rows in x does not match with the KSFFT object, nrows = %d', nrows, n123);
         end
         ng = length(idxnz);
         y = zeros(ng,ncols);
         for j = 1:ncols 
            a3d = fftn(reshape(x(:,j),n1,n2,n3));
            y(:,j) = a3d(idxnz);
         end;
         y = y*vol/n123;
      else
         error('KSFFT: something wrong with the FFT configuration');
      end
   else
      error('KSFFT must be applied to numeric data');
   end;
else
   error('KSFFT syntax: y=F*x')
end;
