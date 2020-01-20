function [s,ext] = kssolvpptype(s,ext)
% KSSOLVPPTYPE Pseudo Potential type for KSSOLV.
%   s = KSSOLVPPTYPE returns a string that is the current desired pseudo
%   potential type in KSSOLV. The default return is string 'default'.
%
%   [s,ext] = KSSOLVPPTYPE returns a string, s, that is the current desired
%   pseudo potential type in KSSOLV together with the extension type, ext.
%   The default extension type is 'upf'.
%
%   KSSOLVPPTYPE(s) sets the pseudopotential file type to be s.
%
%   KSSOLVPPTYPE(s,ext) sets the pseudopotential file type to be s and the
%   extension type to be ext.
%
%   See also matlabroot.
%

%  Copyright (c) 2016-2017 Yingzhou Li and Chao Yang,
%                          Stanford University and Lawrence Berkeley
%                          National Laboratory
%  This file is distributed under the terms of the MIT License.

global kssolvpptype_global kssolvpptypeext_global;

if nargin == 1
    kssolvpptype_global = s;
    kssolvpptypeext_global = 'upf';
elseif nargin == 2
    kssolvpptype_global = s;
    kssolvpptypeext_global = ext;
end

if isempty(kssolvpptype_global)
    kssolvpptype_global = 'default';
end

if isempty(kssolvpptypeext_global)
    kssolvpptypeext_global = 'upf';
end

s = kssolvpptype_global;
ext = kssolvpptypeext_global;

end
