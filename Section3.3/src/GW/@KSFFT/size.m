function [nrow, ncol] = size(F)
%
% usage: [nrow, ncol] = size(F);
%
nrow = 0;
ncol = 0;
if (F.inverse)
   nrow = F.n1*F.n2*F.n3;
   ncol = length(F.idxnz);
elseif (F.forward)
   nrow = length(F.idxnz); 
   ncol = F.n1*F.n2*F.n3;
end
