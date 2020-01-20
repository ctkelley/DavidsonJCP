function [mol, H, X, info, varargout] = relaxatoms(mol, H, X, varargin)
% RELAXATOMS optimize the positions of atoms to minimize total energy. 
%    [mol, H, X, info, varargout] = RELAXATOMS(mol, H, X, mdopt,ksopt) computes
%    a stable structure of atom position. The wholemolecule dynamics progress
%    is done by calculating SCF iteration and minimize total energy
%    by moving atoms alternately. If the molecule is at a unsteable
%    state, the force given by electons will cause nuclear to move, which is
%    treated as lessen energy in Hamiltonian mechanics.
%    The optimization process is done by following method:
%    method = 1: use MATLAB optimization toolbox
%    method = 2: nonlinear conjugate gradient implemented by Amartya
%    method = 3: nonlinear conjugate gradient implemented by Overton
%    method = 4: BFGS package by Overton [part of the 
%                HANSO: Hybrid Algorithm for Non-Smooth Optimization (V2.2)
%                http://www.cs.nyu.edu/overton/software/hanso/]
%    method = 5: Fast Inertial Relaxation Engine (FIRE) 
%                DOI: 10.1103/PhysRevLett.97.170201
%    method = 6: Hybrid of FIRE and methoods 3
%
% Copyright (c) 2018-2019 Bichen Lu and Yingzhou Li,
%                         Fudan University and Duke University

nVarargs = length(varargin);
natoms = sum(mol.natoms);

if isempty(varargin)
    mdopt.int_cord=1;
end
if nVarargs == 2
    if any(strcmp(fieldnames(varargin{1}),'verbose'))
        mdopt=varargin{2};
    else
        mdopt=varargin{1};
    end
else
    if any(strcmp(fieldnames(varargin{1}),'verbose'))
        mdopt.int_cord=1;
    elseif any(strcmp(fieldnames(varargin{1}),'method'))
        mdopt = varargin{1};
    end
end
int_cord=mdopt.int_cord;

% Readjust the positions to fit within the supercell
% NOTE: Assuming cuboidal supercell
xxlim = mol.supercell(1,1);
yylim = mol.supercell(2,2);
zzlim = mol.supercell(3,3);

xtemp = mol.xyzlist(:,1) - (sum(mol.xyzlist(:,1))/natoms - xxlim/2.0); 
ytemp = mol.xyzlist(:,2) - (sum(mol.xyzlist(:,2))/natoms - yylim/2.0); 
ztemp = mol.xyzlist(:,3) - (sum(mol.xyzlist(:,3))/natoms - zzlim/2.0); 

xyztemp = [xtemp ytemp ztemp];

%converte nucleus position into internal-coordinate by shifting the first nucleus to
%position O and rotating the molecule if the int_cord is true
if int_cord == true
    xyz1 = repmat(xyztemp(1,:),natoms,1);
    xyztemp = xyztemp - xyz1;
    zcord=zeros(natoms,3);
    if natoms>=2
        zcord(2,1) = norm(xyztemp(1,:)-xyztemp(2,:));
        zcord(2,2) = 0;
        zcord(2,3) = 0;
        if natoms>=3
            zcord(3,1) = (xyztemp(3,:)*xyztemp(2,:)')/zcord(2,1);
            zcord(3,2) = sqrt(xyztemp(3,:)*xyztemp(3,:)'-zcord(3,1)^2);
            zcord(3,3) = 0;
            if natoms>=4
                zcord(4:end,1) = xyztemp(4:end,:)*xyztemp(2,:)' / zcord(2,1);
                zcord(4:end,2) = (xyztemp(4:end,:)*xyztemp(3,:)'-zcord(3,1)*zcord(4:end,1)) / zcord(3,2);
                zcord(4:end,3) = sqrt(diag(xyztemp(4:end,:)*xyztemp(4:end,:)')-(zcord(4:end,1).^2+zcord(4:end,2).^2));
            end
        end
    end
    xyztemp=zcord;
end
mol = set(mol,'xyzlist',xyztemp);
    
xyzforces = mol.xyzforce;

ksopts = setksopt;
ksopts = setksopt(ksopts,'scftol', 1e-6, 'maxscfiter', 100);
if ( isempty(xyzforces) )
    % run SCF to obtain the forces on the atoms in MOL
    [mol,H,X,~] = scf4m(mol);
%     ksopts = setksopt(ksopts,'rho0',H.rho, 'X0',X, 'verbose','off');      % USE with caution. Seems to worsen SCF convergence 
%     ksopts = setksopt('verbose','off');
    ksopts = setksopt(ksopts,'rho0',H.rho, 'verbose','off');
else
    if (exist('H'))
        ksopts = setksopt(ksopts, 'rho0',H.rho, 'verbose','off');           % USE with caution. Seems to worsen SCF convergence 
%         ksopts = setksopt('verbose','off');
    else
%         ksopts = setksopt(ksopts, 'X0',X, 'verbose','off');       % USE with caution. Seems to worsen SCF convergence 
        ksopts = setksopt('verbose','off');
    end
end

%read mdopt and ksopt or set by default
if isempty(varargin)
    mdopt=setmdopt(mol,ksopts);
end
if nVarargs == 2
    if any(strcmp(fieldnames(varargin{1}),'verbose'))
        ksopts=varargin{1};
    else
        ksopts=varargin{2};
    end
else
    if any(strcmp(fieldnames(varargin{1}),'verbose'))
        ksopts=varargin{1};
        mdopt=setmdopt(mol,ksopts);
    end
end
    
% get the initial force (gradient)

global saved_mol saved_H saved_X;
saved_mol = mol;
saved_H = H;
saved_X = X;

xyzlist = mol.xyzlist;
x0 = xyzlist(:);
if (mdopt.method == 1)
    % requires MATLAB optimization toolbox
    optimopts = mdopt.fminunc_opts;
    
    [x,~] = fminunc(@(x) ksfg(x,mol,ksopts,int_cord), x0, optimopts);
    xyzlist_new = reshape(x, [natoms 3]);
    mol = set(mol,'xyzlist',xyzlist_new);
    if nargout == 6
        [mol,~, ~, info] = scf4m(mol, ksopts,int_cord);
        Etot = info.Etot;
        varargout{1} = Etot;
        varargout{2} = ksfg(x,mol,ksopts,int_cord);
    end
    
elseif (mdopt.method == 2)
    %use nonlinear CG by Amartya
    NLCG_opts = mdopt.NLCG_opts;
    [mol, FORCE, Etot] = NLCG(mol, NLCG_opts, int_cord);
    if nargout == 6
        varargout{1} = Etot;
        varargout{2} = FORCE;
    end
    
elseif (mdopt.method == 3)
    % default: nonlinear CG by M. Overton
    nlcg_pars = mdopt.nlcg_pars;
    nlcg_opts = mdopt.nlcg_opts;
    
    % output order: [x, f, g, frec, grec, alpharec]
    [x, f, g, frec, grec, ~] = nlcg(nlcg_pars, nlcg_opts);

    nIonicIter = length(frec{1}(:));
    natoms = sum(mol.natoms);
    
    if nargout == 8
        varargout{1} = f;                           % converged energy value
        varargout{2} = -reshape(g,[natoms 3]);      % converged force value
        varargout{3} = frec;        % history of energy values
        
        % Get the history of force on each atom at every iteration:
        fAtoms = cell(nIonicIter,1);
        
        for n = 1:nIonicIter
            fAtoms{n} = -reshape(grec{1}(:,n),[natoms 3]);
        end
        varargout{4} = fAtoms;
    end 
    
    fprintf('******************************************* \n');
    fprintf('Number of Ionic Iterations to Converge: %3d \n', nIonicIter);
    fprintf('******************************************* \n');
    
    xyzlist_new = reshape(x,[natoms 3]);
    mol = set(mol,'xyzlist',xyzlist_new);
    
elseif (mdopt.method == 4)
    % The BFGS package by Michael Overton (part of the
    % HANSO: Hybrid Algorithm for Non-Smooth Optimization (V2.2)
    % http://www.cs.nyu.edu/overton/software/hanso/
    %
        % output order: [x, f, d, H, iter, info, X, G, w, fevalrec, xrec, drec, Hrec]
    BFGS_pars = mdopt.BFGS_pars;
    BFGS_opts = mdopt.BFGS_opts;
    BFGS_pars.int_cord=int_cord;
    
    % Call BFGS
    % output order: [x, f, d, H, iter, info, X, G, w, fevalrec, xrec, drec, Hrec]
    [x, f, d, ~, nIonicIter, info, ~, ~, ~, frec, ~, drec, ~] = bfgs(BFGS_pars, BFGS_opts);

    if nargout == 7
        varargout{1} = f;                           % converged energy value
        varargout{2} = reshape(x,[natoms 3]);       % converged position
        varargout{3} = frec;                        % history of energy values
    end
    
    if nargout == 8
        varargout{1} = f;                           % converged energy value
        varargout{2} = -reshape(d,[natoms 3]);      % converged force value
        varargout{3} = frec;        % history of energy values
        
        % Get the history of force on each atom at every iteration:
        fAtoms = cell(nIonicIter,1);
        
        for n = 1:nIonicIter
            fAtoms{n} = -reshape(drec{n}(:,1),[natoms 3]);
        end
        varargout{4} = fAtoms;
    end 
    
    fprintf('EXIT INFO STATUS: %d\n', info);
    
    fprintf('******************************************* \n');
    fprintf('Number of Ionic Iterations to Converge: %3d \n', nIonicIter);
    fprintf('******************************************* \n');
    
    natoms = sum(mol.natoms);
    xyzlist_new = reshape(x,[natoms 3]);
    mol = set(mol,'xyzlist',xyzlist_new);
    
    
elseif (mdopt.method == 5)
    fire_pars=mdopt.fire_pars;
    fire_opts=mdopt.fire_opts;
    [mol, ~, POS, ~, FORCE, Etot, INFO, exitFlag] = quenchbyfire(mol,fire_pars,fire_opts,ksopts,int_cord);
    
    if nargout == 9
        varargout{1} = POS;
        varargout{2} = FORCE;
        varargout{3} = Etot;
        varargout{4} = INFO;
        varargout{5} = exitFlag;
    end
   
    fprintf('******************************* \n');
    fprintf('Number of Ionic Iterations: %3d \n', length(Etot));
    fprintf('******************************* \n');
    
elseif (method == 6)    % hybrid with fire
    
    if nVarargs == 4
        
    [mol, FORCE, Etot] = hybridoptimizer(mol, varargin{1}, varargin{2}, ...
                                         varargin{3}, varargin{4});
    end
    
    if nargout == 6
        varargout{1} = FORCE;
        varargout{2} = Etot;
    end

else
    fprintf('Error: invalid optimization method = %d\n', method);
    fprintf('method must be 1, 2, 3, 4, or 5\n');
end
%
% run SCF again to pass out H and X
ksopts = setksopt();
ksopts.maxscfiter = 100;
ksopts.betamix = 0.1;
ksopts.mixdim = 10;
ksopts.scftol = 1e-6;

ksopts = setksopt();
% ksopts = setksopt(ksopts,'mixtype','broyden','eigmethod','eigs');
ksopts = setksopt(ksopts,'mixtype','pulay','eigmethod','eigs');
% ksopts = setksopt;
% 
% ksopts = setksopt(ksopts, ...
%     'maxscfiter',100,'scftol',1e-10,'cgtol',1e-10,'maxcgiter',100);
[mol,H,X,info]=scf4m(mol,ksopts);

end
%-------------------------------------------------------------------
