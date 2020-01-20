function val = get(F,attr_name)
%
% usage: val = get(F,attr_name);
%        retrieve various attributes of an KSFFT object F.

switch attr_name
  case 'n1'
    val = F.n1;
  case 'n2'
    val = F.n2;
  case 'n3'
    val = F.n3;
  case 'idxnz'
    val = F.idxnz;
  case 'forward'
    val = F.forward;
  case 'inverse'
    val = F.inverse;
  otherwise
    error('invalid attribute');
end;
