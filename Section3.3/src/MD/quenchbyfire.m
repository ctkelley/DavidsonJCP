function [mol, POS, POSRecMat, FORCE, FORCERecMat, Etot, INFO, exitFlag] = quenchbyfire(mol, pars, opts, scfOpts,int_cord)
% 
% Subhajit Banerjee
% June 2017
% SSG@CRD, LBNL
% 

    %fprintf('FIRE executed with default parameter values, options, and mass = 4 \n');
    % 
    % For appropriate mass value please review the following reference
    % http://nanosurf.fzu.cz/wiki/doku.php?id=fire_minimization
    % Also the mass is same for all atoms in FIRE
    %
mass = pars.mass;     % MASS = 4 unit
    %
    % Source:
    % https://journals.aps.org/prl/abstract/10.1103/PhysRevLett.97.170201
    %
    % Qualitatively they should satisfy: nMin larger than 1 (at
    % least a few smooth steps after freezing); fInc larger than but 
    % near to one (avoid too fast acceleration); fDec smaller than 
    % 1 but much larger than zero (avoid too heavy slowing down),
    % alphaStart larger than, but near to zero (avoid too heavy damping); 
    % fAlpha smaller than, but near to one (mixing is efficient some time 
    % after restart).
    %

    %
    % 1 Femtoseconds [fs] = 41.3413745758 atomic time unit
    % dt = 41.3413745758;         % in atomic time unit

    % 
    % Set up options for Kohn-Sham Eq.
    %
if ~exist('scfOpts')
    scfOpts = setksopt();
    scfOpts.maxscfiter = 100;
    scfOpts.betamix = 0.1;
    scfOpts.mixdim = 10;
    scfOpts.scftol = 1.0E-4;
    scfOpts = setksopt(scfOpts,'mixtype','pulay','eigmethod','eigs');
end

% 
dt = opts.dt;
MAXITER = opts.MAXITER;
TOL = opts.TOL;
% 
% Initialization
% 
nAtoms = sum(mol.natoms);
xyzlist = mol.xyzlist;
x0 = xyzlist(:);
v0 = zeros(size(x0));               % velocity initialized to 0 for FIRE
alpha = pars.alphaStart;            % needed for the first FIRE run

% POS = cell(MAXITER, 1);
% FORCE = cell(MAXITER, 1);
% Etot = zeros(MAXITER,1);
% INFO = cell(MAXITER, 1);
if ~exist('int_cord')
    int_cord=true;
end
    
verletParams = struct('v0',v0, 'dt',dt, 'mass',mass);  % initial parameters for Velocity Verlet
% 
% First velocity update with initial velocity = 0
%
fprintf('\n');
fprintf('First velocity Verlet update with starting velocity = 0 \n');
fprintf('\n');
% 
% [x, v, F, mol, H, X, info] = velverletfire(mol, verletParams, scfOpts);
[x, v, F, mol, H, ~, info] = velverletfire(mol, verletParams, scfOpts, int_cord);
% 
it = 0;
cut = 0;
exitFlag = false;

while (it <= MAXITER)
    it = it + 1;
    Etot(it) = info.Etot;
    POS{it} = reshape(x, [nAtoms 3]);
    FORCE{it} = reshape(F, [nAtoms 3]);
    INFO{it} = info;
    
    % Check for convergence (alternative convergence criterion)
    % fCheck = mol.xyzforce;          % [natoms x 3] matrix == reshaped form
    %                                 % of F as output by velverletfire()
    % magF = zeros(size(fCheck,1));
    % for i = 1:length(magF),
    %    magF(i) = norm(fCheck(i,:));
    % end
    % checkConv = max(magF);
    %        
    checkConv = norm(F);
    
    if checkConv <= TOL
        exitFlag = true;
        break
    else
        xyzlist_new = reshape(x,[nAtoms 3]);
        mol = set(mol,'xyzlist',xyzlist_new);
        % scfOpts = setksopt(scfOpts,'rho0',H.rho, 'X0',X, 'verbose','off');% USE with caution. Seems to worsen SCF convergence 
        scfOpts = setksopt(scfOpts,'rho0',H.rho, 'verbose','off');  % USE with caution. Seems to worsen SCF convergence       
%         scfOpts = setksopt(scfOpts,'verbose','off');
    
        % Update v, dt, alpha using FIRE
        fprintf('\n');
        fprintf('Updating v, dt, and alpha using FIRE for iteration: %-5d\n', it);
        fprintf('\n');
        [v, dt, alpha, cut] = fireupdater(it, dt, F, v, cut, alpha, pars);
        % 
        % rerun Velocity Verlet with initial velocity from FIRE update
        % 
        % reassigning parameters for subsequent velocity Verlet run as per
        % the last FIRE output. This changes v0 nd dt.
        verletParams = struct('v0',v, 'dt',dt, 'mass',mass);
        
        fprintf('\n');
        fprintf('Starting velocity Verlet update from FIRE output for iteration: %-5d\n', it);
        fprintf('\n');
        %
        [x, v, F, mol, H, X, info] = velverletfire(mol, verletParams, scfOpts,int_cord);
        if info.converge == 1
            fprintf('\n');
            fprintf('SCF Converged in %d iterations\n', length(info.SCFerrvec));
            fprintf('\n');
        else
            fprintf('\n');
            fprintf('SCF did not converge. Simulation continued to the next (%-5d) FIRE iteration \n', it + 1);
            fprintf('\n');
        end
        
        fprintf('\n');
        fprintf('Completed velocity Verlet update for iteration: %-5d\n', it);
        fprintf('\n');     
        fprintf('\n Norm of force vectro = %-16.10e.\n\n', norm(F));
        % fprintf('\n Maximum force component = %-16.10e.\n\n', max(abs(F)));
        fprintf('\n'); 
        
        % 
        if it == MAXITER
            fprintf('\n');
            fprintf(' FIRE did not converge! \n');
            fprintf('\n');
            fprintf(' Reached maximum number ( %-5d) of iterations! \n', it);
            fprintf('\n');
        end
        
    end
end

if exitFlag == 1
    % maxNnzEl = it;
    % Etot =  Etot(1:maxNnzEl);
    % INFO = INFO{1:maxNnzEl};
    
    fprintf('\n');
    fprintf(' FIRE converged in %-5d steps \n', it);
    fprintf(' The local minumum achieved is = %-16.10f \n', Etot(end));
    fprintf('\n');
    
    % POSRecMat = [POS{1:maxNnzEl}];
    POSRecMat = [POS{:}];       % blocks of [nAtoms 3] arrays arranged 
                                % in a row. total of 'it' blocks, 
                                % one for each iterations
                                
    % POSValAtMin = POSRecMat(:,end-2:end);
    POSValAtMin = [POS{end}];
    
    mol = set(mol,'xyzlist',POSValAtMin);
    
    % FORCERecMat = [FORCE{1:maxNnzEl}];    
    FORCERecMat = [FORCE{:}];    % blocks of [nAtoms 3] arrays arranged 
                                 % in a row. total of 'it' blocks, 
                                 % one for each iterations
                                    
	% FORCEValAtMin = FORCERecMat(:,end-2:end);
    FORCEValAtMin = [FORCE{end}];
    
    fprintf(' The minimum is achieved for Configuration (argmin) \n');
    fprintf('\n');
    disp(POSValAtMin(:,:));
    fprintf('\n');
    fprintf( 'The corresponding forces: \n');
    fprintf('\n');
    % disp(FORCEValAtMin(:,end-2:end)); 
    disp(FORCEValAtMin); 
    fprintf('\n');
    
else
    [minEtot, ind] = min(Etot);
    fprintf('\n');
    fprintf(' The infimum achieved is = %-16.10f \n', minEtot);
    fprintf('\n');
    
    POSRecMat = [POS{:}];       % blocks of [nAtoms 3] arrays arranged 
                                % in a row. total of MAXITER blocks, 
                                % one for each iterations

    POSValAtInf = POSRecMat(:,(ind-1)*3+1:ind*3);
    
    mol = set(mol,'xyzlist',POSValAtInf);
    
    FORCERecMat = [FORCE{:}];    % blocks of [nAtoms 3] arrays arranged 
                                 % in a row. total of MAXITER blocks, 
                                 % one for each iterations
                                    
	FORCEValAtInf = FORCERecMat(:,(ind-1)*3+1:ind*3);
    
    fprintf(' The infimum is achieved for Configuration (argmin) \n');
    fprintf('\n');
    disp(POSValAtInf(:,:));
    fprintf('\n');
    fprintf('at the iteration # %-5d \n', ind);
    fprintf('\n');
    fprintf( 'The corresponding forces: \n');
    fprintf('\n');
    % disp(FORCEValAtInf(:,end-2:end)); 
    disp(FORCEValAtInf);
    fprintf('\n');
end
 
end

function [x, v, F, mol, H, X, info] = velverletfire(mol, params, scfOpts,int_cord)
% 
% Subhajit Banerjee
% June 2017
% SSG@CRD, LBNL
% 

v0 = params.v0;
dt = params.dt;
mass = params.mass;

% Readjust the positions to fit within the supercell
% NOTE: Assuming cuboidal supercell
natoms = sum(mol.natoms);
xxlim = mol.supercell(1,1);
yylim = mol.supercell(2,2);
zzlim = mol.supercell(3,3);

mol_xyzlist = mol.xyzlist;

for i = 1:natoms
    if mol_xyzlist(i,1) < 0.0
        mol_xyzlist(i,1) = mol_xyzlist(i,1) + xxlim;
    elseif mol_xyzlist(i,1) > xxlim
        mol_xyzlist(i,1) = mol_xyzlist(i,1) - xxlim;
    end
    
    if mol_xyzlist(i,2) < 0.0
        mol_xyzlist(i,2) = mol_xyzlist(i,2) + yylim;
    elseif mol_xyzlist(i,2) > yylim
        mol_xyzlist(i,2) = mol_xyzlist(i,2) - yylim;
    end
    
    if mol_xyzlist(i,3) < 0.0
        mol_xyzlist(i,3) = mol_xyzlist(i,3) + zzlim;
    elseif mol_xyzlist(i,3) > zzlim
        mol_xyzlist(i,3) = mol_xyzlist(i,3) - zzlim;
    end
end

mol = set(mol,'xyzlist',mol_xyzlist);

xyzlist = mol.xyzlist;          % [natoms x 3] matrix
x0 = xyzlist(:);                % a vector now
F0 = mol.xyzforce;              % F0 is force (a [nAtoms x 3] matrix)  

F0 = F0(:);                     % F0 is a vector now
%
% Position update in Velocity Verlet
% 
x = x0 + v0*dt + 0.5*dt*dt*F0/mass;              % Note: F0 is acceleration

% Update the position of molecule object
nAtoms = sum(mol.natoms);
xyzlist_new = reshape(x,[nAtoms 3]);
mol = set(mol,'xyzlist',xyzlist_new);

[mol, H, X, info] = scf4m(mol, scfOpts, int_cord);     % This should be the only SCF cycle
 
F = mol.xyzforce;               % F is force (a [nAtoms x 3] matrix) 

F = F(:);                       % F is a vector now

% Velocity update in Velocity Verlet
v = v0 + 0.5*(F0 + F)/mass*dt;

end