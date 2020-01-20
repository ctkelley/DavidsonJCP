function [f,g] = ksef(x,pars)
global saved_mol saved_H saved_X
mol = pars.mol;
if ~exist('pars.int_cord')
    int_cord=true;
else
    int_cord=pars.int_cord;
end
ksopts = pars.ksopts;
natoms = sum(mol.natoms);
xyz = reshape(x, [natoms 3]);
mol = set(mol,'xyzlist',xyz);
if (~isempty(saved_H) && ~isempty(saved_X))
    % ksopts = setksopt(ksopts,'X0',saved_X,'rho0', saved_H.rho);   % USE with caution. Seems to worsen SCF convergence
    ksopts = setksopt(ksopts,'rho0', saved_H.rho);
end
[saved_mol, saved_H, saved_X, info] = scf4m(mol, ksopts,int_cord);
f = info.Etotvec(end);
xyzforces = saved_mol.xyzforce;
g = -xyzforces(:); % no negative sign here
fprintf('norm(g) = %11.3e\n', norm(g));
end