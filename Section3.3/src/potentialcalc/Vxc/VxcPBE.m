function [vxc,uxc,rho] = VxcPBE(mol,rho)
% VXCPBX PBE exchange correlation.
%    [vxc,uxc,rho] = VXCPBE(mol,rho) returns the PBE exchange correlation
%    of the rho.
%
%   See also exRef.

%  Copyright (c) 2016-2017 Yingzhou Li and Chao Yang,
%                          Stanford University and Lawrence Berkeley
%                          National Laboratory
%  This file is distributed under the terms of the MIT License.

e2 = e2Def();

[grho,grho2] = getgradrho(mol,rho);

vxc = zeros(size(rho));
uxc = zeros(size(rho));

idxxc = abs(rho) > 1e-10;
rhoxc = abs(rho(idxxc));

rs = ((3/(4*pi))./rhoxc).^(1/3);
[vx,ux] = VExchange_sla(rs);
[vc,uc] = VCorrelation_pw(rs);

vxc(idxxc) = e2*(vx+vc);
uxc(idxxc) = e2*(ux+uc);

idxcxc = abs(rho) > 1e-6 & grho2 > 1e-10;
rhocxc = abs(rho(idxcxc));
grho2cxc = grho2(idxcxc);

[v1gcx,v2gcx,ugcx] = VGCExchange_pbx(rhocxc,grho2cxc);
[v1gcc,v2gcc,ugcc] = VGCCorrelation_pbc(rhocxc,grho2cxc);

vxc(idxcxc) = vxc(idxcxc) + e2*(v1gcx+v1gcc);
uxc(idxcxc) = uxc(idxcxc) + e2*(ugcx+ugcc);

h2xc = zeros(size(rho));
h2xc(idxcxc) = e2*(v2gcx+v2gcc);
h2xcgrho2 = repmat((h2xc),1,1,1,3).*grho;   % this would issue warning 
                                            % for all the versions of 
                                            % MATLAB below R2013b. See
                                            % http://www.dynare.org/DynareWiki/MatlabVersionsCompatibility
                                            % for details.
n1 = mol.n1;
n2 = mol.n2;
n3 = mol.n3;
ecut2 = mol.ecut2;
grid2 = Ggrid(mol,ecut2);
idxnz2 = grid2.idxnz;
gkx2 = grid2.gkx;
gky2 = grid2.gky;
gkz2 = grid2.gkz;
hxc2g = reshape(fft3(h2xcgrho2),n1*n2*n3,3);
gaux = 1i*(hxc2g(idxnz2,1).*gkx2 ...
    + hxc2g(idxnz2,2).*gky2 + hxc2g(idxnz2,3).*gkz2);
gauxfull = zeros(n1*n2*n3,1);
gauxfull(idxnz2) = gaux;
v2xc = real(ifft3(reshape(gauxfull,n1,n2,n3)));

vxc = vxc-v2xc;

end
