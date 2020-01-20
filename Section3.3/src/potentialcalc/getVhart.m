function vhart = getVhart(mol,rho)

% extracting the discretiation dimension in real space
n1 = mol.n1;
n2 = mol.n2;
n3 = mol.n3;
% energy cut-off from the molecule
ecut2 = mol.ecut2;
% add temperarily
grid2  = Ggrid(mol,ecut2);
% indeces of the non-zero elements in Fourier space
idxnz2 = grid2.idxnz;
% Fourier grid
gkk2   = grid2.gkk;

% Calculate Hartree potential
rhog   = fft3(rho);
rhog   = rhog(idxnz2);

w      = zeros(n1,n2,n3);
inz    = abs(gkk2) ~= 0;
w(idxnz2(inz)) = e2Def()*4*pi*rhog(inz)./gkk2(inz);
vhart  = real(ifft3(w));

end
