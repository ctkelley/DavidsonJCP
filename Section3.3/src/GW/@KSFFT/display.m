function display(F)
%
% display method for an KSFFT object
%
fprintf('KSFFT configuration: \n');
fprintf('n1 = %d, n2 = %d, n3 = %d\n', F.n1, F.n2, F.n3);
fprintf('ng = %d\n', length(F.idxnz));
if (F.forward) 
   fprintf('forward\n');
else
   fprintf('inverse\n');
end;
