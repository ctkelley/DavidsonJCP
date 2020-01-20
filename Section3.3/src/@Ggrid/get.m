function val = get(a,attr_name)
%
% usage: val = get(a,attr_name);
%        retrieve various attributes of an FreqMask object a.
%
switch attr_name
  case 'ecut'
    val = a.ecut;
  case 'ng'
    val = a.ng;
  case 'idxnz'
    val = a.idxnz;
  case 'gkk'
    val = a.gkk;
  case 'gkx'
    val = a.gkx;
  case 'gky'
    val = a.gky;
  case 'gkz'
    val = a.gkz;
  otherwise
    error('Invalid attribute');
end;
