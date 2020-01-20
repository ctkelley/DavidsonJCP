function F = set(F,varargin)
%
% usage: F = set(F, attr_name1, value1, attr_name2, value2,...);
% set the attributes of a KSFFT object F.
%
Hargin = varargin;
while (length(Hargin)>=2),
   attr_name = Hargin{1};
   value     = Hargin{2};
   Hargin    = Hargin(3:end);
   switch attr_name
     case 'n1'
        F.n1 = value;
     case 'n2'
        F.n2 = value; 
     case 'n3'           
        F.n3 = value;    
     case 'idxnz'
        F.idxnz = value;
     case 'forward'
        F.forward= value;
     case 'inverse'
        F.inverse = value;
        % 
     otherwise
        fprintf('Attribute name: %s\n', attr_name);
        error('cannot be found');
   end;
end;
