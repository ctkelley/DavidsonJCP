function val = get(a,attr_name)
%
% usage: val = get(a,attr_name);
%        retrieve various attributes of an Atom object a.
% e.g.  asymbol = get(a,'symbol') returns the atomic symbol associated with 
%       the Atom object a.
%       an = get(a,'anum') returns the atomic number associated with a
switch attr_name
  case {'name','Name'}
    val = a.name;
  case 'anum'
    val = a.anum;
  case 'amass'
    val = a.amass;
  case 'venum'
    val = a.venum;
  otherwise 
    error('invalid attribute requested');
end;

