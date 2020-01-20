function [vlocg, rhog] = vloc2g(pp,gl)
% VLOC2G converts Vloc at irregular grid in real-space to uniform grid in
%        G-space.
%    [vlocg,rhog] = VLOC2G(pp,gl) calculate the Fourier transform of the
%    Vloc in the pseudo potential to G-space at location gl. The input of
%    this function is the pseudo potential and gl. And the function returns
%    the Vloc at G-space on gl and the corresponding rho. The detailed
%    implementation of the Fourier transform are detailed in the documents.
%
%    See also VNL2G, GETVION. 

%  Copyright (c) 2015-2016 Yingzhou Li and Chao Yang,
%                          Stanford University and Lawrence Berkeley
%                          National Laboratory
%  This file is distributed under the terms of the MIT License.

e2 = e2Def();

n = length(pp.r);
ngl = length(gl);
r = pp.r;
rab = pp.rab;
vloc = pp.vloc;
rho = pp.rhoatom;
zp = pp.venum;


%preallocation
vlocg = zeros(ngl,1);
rhog = zeros(ngl,1);

% calculate vlocg
if gl(1) < 1e-8
    aux = r.*(r.*vloc + zp*e2);
    vlocg(1) = simpson(n,aux,rab);
    igl0 = 2;
else
    igl0 = 1;
end

aux1 = r.*vloc+zp*e2*erf(r);
fac = zp*e2;

igl = igl0:ngl;
len = length(igl);
gx = sqrt(gl(igl));
aux = repmat(aux1',len,1).*sin(gx*r')./repmat(gx,1,n);
vlcp = simpson(n,aux,rab,2);
vlocg(igl) = vlcp - fac*exp(-gl(igl)/4)./gl(igl);

% calculate rhog
if gl(1) < 1e-8
    aux = rho;
    rhog(1) = simpson(n,aux,rab);
    igl0 = 2;
else
    igl0 = 1;
end

igl = igl0:ngl;
len = length(igl);
gx = sqrt(gl(igl));
aux = repmat(rho',len,1).*sin(gx*r')./(gx*r');
rzeroidx = r<1e-8;
aux(:,rzeroidx) = repmat(rho(rzeroidx)',len,1);
rhog(igl) = simpson(n,aux,rab,2);

end
