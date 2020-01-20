function [f, g] = ksfg(x, mol, ksopts,int_cord)
global saved_mol saved_H saved_X;
natoms = sum(mol.natoms);
xyz = reshape(x, [natoms 3]);
mol = set(mol,'xyzlist',xyz);
if ( ~isempty(saved_X) && ~isempty(saved_H) )
%     ksopts = setksopt(ksopts,'X0',saved_X,'rho0', saved_H.rho);   % USE with caution. Seems to worsen SCF convergence 
    ksopts = setksopt(ksopts,'rho0', saved_H.rho);
%     disp('use saved')
%     ksopts = setksopt('verbose','off');
end
[saved_mol,saved_H,saved_X,info] = scf4m(mol, ksopts,int_cord);
f = info.Etotvec(end);
xyzforces = saved_mol.xyzforce;
g = -xyzforces(:);
fprintf('norm(g) = %11.3e\n', norm(g));
end