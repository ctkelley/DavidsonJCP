function [cry,BH,BX,info] = scf4c(cry,options)
% SCF4C Self Consistent Field iteration for crystal.
%    [cry,BH,BX,info] = SCF4C(cry,options) adopts Self Consistent Field
%    (SCF) iteration to find the ground state minimum total energy and the
%    corresponding wave functions. cry is a Crystal object and options is
%    the options for running the SCF4C. Please read setksopt for detailed
%    information about options. SCF4C returns the crystal cry, the Bloch
%    Hamiltonian BH, the Bloch wave functions BX, and the information about
%    the iteration. We also support non-zero temperature for the
%    calculation.
%
%   See also scf4m, scf4c, dcm, trdcm.

%  Copyright (c) 2016-2017 Yingzhou Li and Chao Yang,
%                          Stanford University and Lawrence Berkeley
%                          National Laboratory
%  This file is distributed under the terms of the MIT License.

if (nargin < 2)
    options = setksopt();
end

% Set timer
tstart  = cputime;

% Initialize input variables
verbose    = ~strcmpi(options.verbose,'off');
force      = options.force;
maxscfiter = options.maxscfiter;
scftol     = options.scftol;
what2mix   = options.what2mix;
mixtype    = options.mixtype;
mixdim     = options.mixdim;
betamix    = options.betamix;
brank      = options.brank;
BX         = options.X0;
rho        = options.rho0;

nkpts      = cry.nkpts;
wks        = cry.wks;
nspin      = cry.nspin;
temperature= cry.temperature;
Tbeta      = temperature*8.6173324e-5/13.6;

% Initialize Hamiltonian, Wavefun, and Preconditioners
[cry,BH,BX,BHprec,noccs] = iterinit(cry,rho,BX);

% calculate Ewald and Ealphat
Eewald     = getEewald(cry);
Ealphat    = getEalphat(cry);

vion       = BH.vion;
vext       = BH.vext;
vtot       = BH.vtot;
rho        = BH.rho;
nBXcols    = ncols(BX);

% Initialize output variables
Etotvec    = zeros(maxscfiter,1);
scferr     = zeros(maxscfiter,1);
ev         = zeros(sumel(nBXcols),1);
dfmat      = [];
dvmat      = [];
cdfmat     = [];

fprintf('Beging SCF4C calculation for %s...\n',cry.name);
for iterscf = 1:maxscfiter
    
    fprintf('SCF iter %3d:\n', iterscf);
    
    rhoin  = rho;
    vtotin = vtot;
    
    idx = 0;
    for ik = 1:nkpts
        idx = idx(end) + (1:nBXcols(ik));
        [X, ev(idx)] = updateX(cry, BH{ik}, BX{ik}, BHprec{ik}, options);
        BX{ik} = X;
    end
    info.Eigvals = ev;
    
    % TODO: which one should appear first, getocc or multiply with wks?
    %[idxev,ev] = sortev(ev,cry.nkpts*cry.nel/2*cry.nspin);
    [occs,efermi] = getocc(ev,noccs,Tbeta);
    
    idx = 0;
    for ik = 1:nkpts
        idx = idx(end) + (1:nBXcols(ik));
        ev(idx) = ev(idx)*wks(ik);
    end
    
    % Update density function rho
    rho = getcharge(cry,BX,occs);
    BH.rho = rho;
    if strcmpi(what2mix,'rho')
        rhoerr = norm(rho(:)-rhoin(:))/norm(rhoin(:));
        scferr(iterscf) = rhoerr;
        fprintf('Rel Rho Err     = %20.3e\n',rhoerr);
        [rho,dfmat,dvmat,cdfmat] = potmixing(cry,rhoin,rho,...
            iterscf,mixtype, betamix, ...
            dfmat, dvmat, cdfmat, mixdim, ...
            brank);
    end
    
    Entropy = getEntropy(occs,Tbeta);
    
    % Kinetic energy and some additional energy terms
    Ekin = (2/nspin)*sum(ev.*occs);
    
    % ionic and external potential energy was included in Ekin
    % along with incorrect Ecoul and Exc. Need to correct them
    % later;
    Ecor = getEcor(cry, rho, vtot, vion, vext);
    
    % Compute Hartree and exchange correlation energy and potential
    % using the new charge density; update the total potential
    [vhart,vxc,uxc2,rho]=getVhxc(cry,rho);
    
    % Update total potential
    vtot = getVtot(cry, vion, vext, vhart, vxc);
    if strcmpi(what2mix,'pot')
        vtoterr = norm(vtot(:)-vtotin(:))/norm(vtotin(:));
        scferr(iterscf) = vtoterr;
        fprintf('Rel Vtot Err    = %20.3e\n',vtoterr);
        [vtot,dfmat,dvmat,cdfmat] = ...
            potmixing(cry,vtotin,vtot,iterscf,mixtype,...
            betamix,dfmat,dvmat,cdfmat,mixdim,brank);
    end
    BH.vtot = vtot;
    
    % Calculate the potential energy based on the new potential
    Ecoul = getEcoul(cry,abs(rho),vhart);
    Exc   = getExc(cry,abs(rho),uxc2);
    Etot  = Entropy + Ekin + Eewald + Ealphat + Ecor + Ecoul + Exc;
    Etotvec(iterscf) = Etot;
    
    % Convergence check
    fprintf('Total Energy    = %20.13e\n', Etot);
    [cvg,resfro] = reportconverge(BH,BX,iterscf,maxscfiter, ...
        scferr(iterscf),scftol,verbose);
    if cvg
        info.converge = true;
        break;
    end
    
end

BX = assignoccs(BX,occs);
if force
    cry.xyzforce = getFtot(cry,BH,BX,rho);
end

info.Etotvec = Etotvec(1:iterscf);
info.SCFerrvec = scferr(1:iterscf);
info.efermi = efermi;

timetot = cputime - tstart;
fprintf('Etot            = %20.13e\n', Etot);
fprintf('Entropy         = %20.13e\n', Entropy);
fprintf('Ekin            = %20.13e\n', Ekin);
fprintf('Eewald          = %20.13e\n', Eewald);
fprintf('Ealphat         = %20.13e\n', Ealphat);
fprintf('Ecor            = %20.13e\n', Ecor);
fprintf('Ecoul           = %20.13e\n', Ecoul);
fprintf('Exc             = %20.13e\n', Exc);
fprintf('Efermi          = %20.13e\n', efermi);
fprintf('--------------------------------------\n');
fprintf('Total time used = %20.3e\n', timetot);
fprintf('||HX-XD||_F     = %20.3e\n', resfro);

end
