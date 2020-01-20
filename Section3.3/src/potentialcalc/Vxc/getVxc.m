function [vxc,uxc,rho] = getVxc(mol,rho)
% GETVXC exchange correlation.
%    [vxc,uxc,rho] = GETVXC(mol,rho) returns the exchange correlation of
%    the rho. The type of the exchange correlation is determined by the
%    pseudopotential.
%
%   See also exRef.

%  Copyright (c) 2016-2017 Yingzhou Li and Chao Yang,
%                          Stanford University and Lawrence Berkeley
%                          National Laboratory
%  This file is distributed under the terms of the MIT License.

switch upper(mol.ppvar.funct)
    case 'PBE'
        [vxc,uxc,rho] = VxcPBE(mol,rho);
    case 'SLA-PW-PBX-PBC'
        [vxc,uxc,rho] = VxcPBE(mol,rho);
    case 'PZ'
        [vxc,uxc,rho] = VxcPZ(rho);
    otherwise
        [vxc,uxc,rho] = VxcPZ(rho);
end

end