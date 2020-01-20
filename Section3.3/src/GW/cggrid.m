function ga = cGgrid(Ggrid,mol,idx)
%
% ga = cGgrid(Ggrid,mol,idx)
%
% display all or a subset of reciprocal space grid points within the 
% Ecut limit
%
% prints out igkx, igky, igkz, igkk
%
C = get(mol,'supercell');
Lx = norm(C(:,1));
Ly = norm(C(:,2));
Lz = norm(C(:,3));

igkx = Ggrid.gkx*Lx/(2*pi);
igky = Ggrid.gky*Ly/(2*pi);
igkz = Ggrid.gkz*Lz/(2*pi);
igkk = igkx.^2 + igky.^2 + igkz.^2;

if (nargin < 3)
   ga = [igkx igky igkz igkk]
else
   ga = [igkx(idx) igky(idx) igkz(idx) igkk(idx)];
end
