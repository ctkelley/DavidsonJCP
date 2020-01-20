function vnew = kerkmix(mol, vin, vout)
%
% Usage: vnew = kerkmix(mol, vin, vout);
%
% Purpose: 
%     Construct a new 3D potential by mixing the potentials vin and vout
%     using Kerker's recipe;
%
% Input:
%     mol  - a Molecule object
%     vin  - a 3D potential
%     vout - a 3D potential
%
% Output:
%     vnew - a 3D output potential  
%
[n1,n2,n3]=size(vin);
n123 = n1*n2*n3;
vres = vout - vin;
rfft = fftn(vres);
ecut2 = mol.ecut2;
grid2 = Ggrid(mol,ecut2);
inz  = grid2.idxnz;
gkk2 = grid2.gkk;
%
rfft(inz) = (0.8*(gkk2/2)./((gkk2/2)+0.5)-0.8).*rfft(inz);
%
izr  = find(gkk2==0); 
rfft(inz(izr)) = -rfft(inz(izr))/2.0;
%
vcor = ifftn(rfft);
vnew = vin + vcor + 0.8*vres;
