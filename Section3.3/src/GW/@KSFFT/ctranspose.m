function Y = ctranspose(X)
%
% Use to indicate the Fourier transform is a forward transform
% i.e., from real space to reciprocal space
%
% usage: Y = X';
%
if (nargin == 1)
   if (X.forward) 
      X.forward = 0;
      X.inverse = 1;
   elseif (X.inverse)
      X.inverse = 0;
      X.forward = 1; 
   else
      error('not sure which direction the original transform goes');
   end;
   Y = X;
else
   error('KSFFT: invalid syntax')
end;
