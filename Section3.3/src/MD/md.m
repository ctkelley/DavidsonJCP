function [POS, VEL, F, eTot, eNKin, INFO] = md(mol, ksopts, H, X, nSteps)
% 
% Subhajit Banerjee
% June 2017
% SSG@CRD, LBNL

if ( isempty(mol.xyzforce) )
    % run initial SCF to obtain the forces on the atoms, etc. from MOL
    [mol, H, X,~] = scf4m(mol, ksopts);    
    ksopts = setksopt(ksopts,'rho0',H.rho, 'X0',X, 'verbose','off');
else    
    if (nargin > 2)
        ksopts = setksopt(ksopts, 'rho0',H.rho, 'verbose','off');
    end
    if (nargin > 3)
        ksopts = setksopt(ksopts,'rho0',H.rho, 'X0',X, 'verbose','off');
    end 
end

nAtoms = sum(mol.natoms);
xyzlist = mol.xyzlist;
x0 = xyzlist(:);

% Extract the masses (in au):
alist = mol.alist;
uniqMass = [mol.atoms.amass];
mass = uniqMass(alist);     % mass of the atoms arranged as per the 
                            % order appearing in atomlist/xyzlist
                            
% Can use Maxwell-Boltzmann velocity distribution (would need k_B and T)
v0 = rand(size(x0));        % randomly initialized velocity

% 1 Femtoseconds [fs] = 41.3413745758 atomic time unit
dt = 41.3413745758;         % in atomic time unit

params = struct('v0',v0, 'dt',dt);      % initial parameters

POS = cell(nSteps, 1);                % DON'T DO THIS
VEL = cell(nSteps, 1);                % DON'T DO THIS
Force = cell(nSteps, 1);            % DON'T DO THIS
INFO = cell(nSteps, 1);             % DON'T DO THIS

eTot = zeros(nSteps,1);
eNKin = zeros(nSteps,1);

for i = 1:nSteps,
    % The only SCF cycle
    [x, v, F, mol, H, X, info] = velverletmd(mol, mass, ...
                                             params, ksopts);
                            
    vMat = reshape(v, [nAtoms 3]);
    % Store values
    POS{i} = reshape(x, [nAtoms 3]);
    VEL{i} = vMat;
    Force{i} = F;
    INFO{i} = info;
    
    % Calculate the nuclear kinetic energy under B-O Approximation
    ke = 0.0;
    for n = 1:nAtoms,
        vN = vMat(n,:);
        ke = ke + mass(n)*dot(vN, vN);
    end

    eTot(i) = info.Etot + 0.5*ke;
    eNKin(i) = ke;
    
    % Update v, H, and X for the next SCF cycle 
    params = struct('v0',v, 'dt',dt);       % set v <-- v0
    ksopts = setksopt(ksopts,'rho0',H.rho, 'X0',X, 'verbose','off');
    
    fprintf('Time = %-7.2f fs (= %-12.8e au) completed \n', double(i), i*dt);
    
end
