function [x, v, F, mol, H, X, info] = velverletmd(mol, mass, params, ksopts)

% Subhajit Banerjee
% June 2017
% SSG@CRD, LBNL

if ( isempty(mol.xyzforce) )
    % run initial SCF to obtain the forces on the atoms, etc. from MOL
    [mol, H, X,~] = scf4m(mol, ksopts);    
    ksopts = setksopt(ksopts,'rho0',H.rho, 'X0',X, 'verbose','off');
end

v0 = params.v0;
dt = params.dt;

xyzlist = mol.xyzlist;          % [natoms x 3] matrix
x0 = xyzlist(:);                % a vector now
F0 = mol.xyzforce;              % [natoms x 3] matrix 

F0 = computeaccln(F0, mass);    % acceleration F0 is a [natoms x 3] matrix
F0 = F0(:);                     % a vector now
%
% Position update in Velocity Verlet
% 
x = x0 + v0*dt + 0.5*dt*dt*F0;              % F0 is accln now

nAtoms = sum(mol.natoms);
xyzlist_new = reshape(x,[nAtoms 3]);
mol = set(mol,'xyzlist',xyzlist_new);

[mol, H, X, info] = scf4m(mol, ksopts);     % The only SCF cycle
 
F = mol.xyzforce;               % F is force (a [nAtoms x 3] matrix) 
F = computeaccln(F, mass);      % F is acceleration now (a [nAtoms x 3] matrix) 
F = F(:);

% Velocity update in Velocity Verlet
v = v0 + 0.5*(F0 + F)*dt;

end


function accln = computeaccln(F, m)

nAtoms = size(F,1);              % number of atoms
accln = zeros(size(F));

for i = 1:nAtoms
    accln(i,:) = F(i,:)/m(i);
end

end