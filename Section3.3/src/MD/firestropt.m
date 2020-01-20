function[x, fRec] = firestropt(mol, scfOpts, pars, opts)
% Subhajit Banerjee
% June 2017
% SSG@CRD, LBNL

if (nargin == 1), 
    %
    % Source:
    % https://journals.aps.org/prl/abstract/10.1103/PhysRevLett.97.170201
    %
    % Qualitatively they should satisfy: Nmin larger than 1 (at
    % least a few smooth steps after freezing); finc larger than but 
    % near to one (avoid too fast acceleration); fdec smaller than 
    % 1 but much larger than zero (avoid too heavy slowing down),  start 
    % larger than, but near to zero (avoid too heavy damping); f  
    % smaller than, but near to one (mixing is efficient some time 
    % after restart).
    % 
    opts.fDec = 0.5;
    opts.fInc = 1.1;
    opts.nMin = 5;
    opts.alphaStart = 0.1;
    opts.fAlpha = 0.99;
    opts.alpha = opts.alphaStart;
    %
    % 1 Femtoseconds [fs] = 41.3413745758 atomic time unit
    % dt = 41.3413745758;         % in atomic time unit
    pars.dt = 41.3413745758;
    pars.MAXITER = 1000;
    pars.TOL = 1.0E-06;
    % 
    % Set up options for Kohn-Sham Eq.
    %
    scfOpts = setksopt();
    scfOpts.maxscfiter = 200;
    scfOpts.betamix = 0.08;
    scfOpts.mixdim = 6;
    scfOpts.scftol = 1.0E-6;
    scfOpts = setksopt(scfOpts,'mixtype','broyden','eigmethod','eigs');
end
%     
fDec = opts.fDec;
fInc = opts.fInc;
nMin = opts.nMin;
alphaStart = opt.alphaStart;
fAlpha = opts.fAlpha;
alpha = opts.alpha;
% 
dt = pars.dt;
MAXITER = pars.MAXITER;
TOL = pars.TOL;
% 
dtMax = 10.0*dt;
% 
% Initialization
% 
nAtoms = sum(mol.natoms);
xyzlist = mol.xyzlist;
x0 = xyzlist(:);
v0 = zeros(size(x0));               % velocity initialized to 0 for FIRE
    
params = struct('v0',v0, 'dt',dt);  % initial parameters for Velocity Verlet
% 
% First velocity update with initial velocity = 0
% 
[x, v, F, mol] = VelocityVerletFIRE(mol, params, scfOpts);
% 
pNegFlag = 0;
ITER = 1;

while (ITER <= MAXITER),
    
    % Check for convergence
    fCheck = mol.xyzforce;          % [natoms x 3] matrix 
    magF = zeros(size(fCheck,1));

    parfor i = 1:length(magF),
        magF(i) = norm(fCheck(i,:));
    end

    checkConv = max(magF);

    if checkConv < TOL,
        fprintf('Converged in %-4d steps \n', ITER);
        return
    end
    
    xyzlist_new = reshape(x,[nAtoms 3]);
    mol = set(mol,'xyzlist',xyzlist_new);
    
    P = dot(F, v);      % Power
    hatF = F/norm(F);   % Fhat (the unit vector)
    %
    % FIRE velocity update formula
    %
    v = (1.0 - alpha)*v + alpha*hatF*norm(v);

    if P <= 0.0,        
        pNegFlag = pNegFlag + 1;
        v = 0;                      % Reset velosity to 0
        alpha = alphaStart;         % reset alpha to alpha_start
        dt = dt*fDec;               % decrease dt
    elseif (P > 0.0) && (pNegFlag > nMin),
        pNegFlag = 0;               % restart the counter to count neg. power        
        dt = min(dt*fInc, dtMax);   % Update dt
        alpha = fAlpha*alpha;       % update alpha
    end
    
    params = struct('v0',v, 'dt',dt);      % reset parameters for MD

    % Rerun MD
    %
    [x, v, F, mol] = VelocityVerletFIRE(mol, params, scfOpts);
    
    ITER = ITER + 1;
end  

if ITER == MAXITER
    error('Reached maximum number of FIRE iterations!');
end

end