function [mol,H,X,info] = trdcm(mol,options)
% TRDCM Trust Region Direct Constrained Minimization algorithm.
%    [mol,H,X,info] = TRDCM(mol,options) adopts Trust Region Direct
%    Constrained Minimization algorithm to find the ground state minimum
%    total energy and the corresponding wave functions. mol is a Molecule
%    object and options is the options for running the TRDCM. Please read
%    setksopt for detailed information about options. TRDCM returns the
%    molecule mol with/without force, the Hamiltonian H, the wave functions
%    X, and the information for each iteration. If a crystal object is
%    input instead of molecule, the results will be generated by trdcm4c
%    function.
%
%   See also scf, dcm, trdcm4m, trdcm4c.

%  Copyright (c) 2016-2017 Yingzhou Li and Chao Yang,
%                          Stanford University and Lawrence Berkeley
%                          National Laboratory
%  This file is distributed under the terms of the MIT License.

if (nargin < 2)
    options = setksopt();
end

if mol.temperature > 0
    [mol,H,X,info] = trdcm4m(mol,options);
    return;
elseif isa(mol,'Crystal')
    [mol,H,X,info] = trdcm4c(mol,options);
    return;
end

% Set timer
tstart  = cputime;

% Initialize input variables
verbose    = ~strcmpi(options.verbose,'off');
force      = options.force;
maxdcmiter = options.maxdcmiter;
maxinerscf = options.maxinerscf;
trdcmtol   = options.dcmtol;
X          = options.X0;
rho        = options.rho0;

maxtry     = 30;
fudge      = trdcmtol;

nspin      = mol.nspin;

% Initialize Hamiltonian, Wavefun, and Preconditioners
[mol,H,X,Hprec,nocc] = iterinit(mol,rho,X);

% calculate Ewald and Ealphat
Eewald     = getEewald(mol);
Ealphat    = getEalphat(mol);

vion       = H.vion;
vext       = H.vext;
vtot       = H.vtot;

% Initialize output variables
Etotvec    = zeros(maxdcmiter,1);
trdcmerr     = zeros(maxdcmiter,1);

fprintf('Beging TRDCM calculation for %s...\n',mol.name);

% One single iteration of scf
[X,ev] = updateX(mol, H, X, Hprec, options);
rho = getcharge(mol,X);
Ekin = (2/nspin)*sum(ev(1:nocc));
Ecor = getEcor(mol, rho, vtot, vion, vext);
[vhart,vxc,uxc2,rho]=getVhxc(mol,rho);
vtot = getVtot(mol, vion, vext, vhart, vxc);
H.vtot = vtot;
Ecoul = getEcoul(mol,abs(rho),vhart);
Exc   = getExc(mol,abs(rho),uxc2);
Etot = Eewald + Ealphat + Ekin + Ecor + Ecoul + Exc;

for iterdcm = 1:maxdcmiter
    fprintf('DCM iter = %d\n', iterdcm);
    
    vtotin = vtot;
    rhoin = rho;
    
    HX = H*X;
    T = X'*HX;
    T = (T+T')/2;
    R = HX - X*T;
    for j = 1:nocc
        R(:,j)=Hprec.*R(:,j);
    end
    
    % construct the subspace
    if iterdcm == 1
        Y = [X R];
    else
        Y = [X R P];
    end
    
    % project the kinetic and ionic potential part of the Hamiltonian
    % into the subspace spanned by Y
    T = Y'*applyKIEP(H,Y);
    
    B = Y'*Y;
    B = (B+B')/2;
    G = eye(size(B));
    G = G(:,1:nocc);
    
    % solve the projected problem
    sigma = 0;
    numtry  = 0;
    for iterscf = 1:maxinerscf
        Etot0   = Etot;
        
        % project the nonlinear potential part of the Hamiltonian
        A = T + Y'*applyNP(H,Y);
        A = (A+A')/2;
        BG = B*G;
        C = BG*BG';
        C = (C+C')/2;
        [G,D]=eig(A-sigma*C,B,'chol');
        ev = diag(D);
        G = G(:,1:nocc);
        
        % update wavefunction and charge density
        X = Y*G;
        rho = getcharge(mol,X);
        
        %  update the the Hartree and exchange correlation potential
        %  based on the new charge density rho
        [vhart,vxc,uxc2,rho]=getVhxc(mol,rho);
        
        % kinetic energy
        Ekin = (2/nspin)*trace( G'*T*G );
        
        % Calculate the potential energy based on the new potential
        Ecoul = getEcoul(mol,abs(rho),vhart);
        Exc   = getExc(mol,abs(rho),uxc2);
        Etot = Eewald + Ealphat + Ekin + Ecoul + Exc;
        %
        if (Etot > Etot0 )
            %
            % total energy increased, impose trust region by adjusting sigma
            %
            gaps   = ev(2:end)-ev(1:end-1);
            gapmax = max(gaps);
            gap0   = ev(nocc+1)-ev(nocc);
            while gap0 < 0.9*gapmax && numtry < maxtry
                sigma = (sigma==0)*2*gapmax + 2*sigma;
                [G,D] = eig(A-sigma*C,B,'chol');
                ev = diag(D);
                gapmax = max(ev(2:end)-ev(1:end-1));
                gap0   = ev(nocc+1)-ev(nocc);
                numtry = numtry + 1;
            end
            G = G(:,1:nocc);
        end
        
        % --- check Etot again ---
        while Etot > Etot0 && abs(Etot-Etot0)>fudge*abs(Etot0) ...
                && numtry < maxtry
            X = Y*G;
            rho = getcharge(mol,X);
            
            [vhart,vxc,uxc2,rho]=getVhxc(mol,rho);
            
            % kinetic energy
            Ekin = (2/nspin)*trace(G'*T*G);
            
            % the potential energy based on the new potential
            Ecoul = getEcoul(mol,abs(rho),vhart);
            Exc   = getExc(mol,abs(rho),uxc2);
            Etot = Eewald + Ealphat + Ekin + Ecoul + Exc;
            
            % increase sigma again if Etot is still larger
            if (Etot > Etot0)
                sigma = (sigma==0)*2*gapmax + 2*sigma;
                BG = B*G;
                C = BG*BG';
                C = (C+C')/2;
                [G,~] = eig(A-sigma*C,B,'chol');
                G = G(:,1:nocc);
                numtry = numtry + 1;
            end
        end
        
        % Update the nonlinear potential only
        H.vnp = vhart+vxc;
        vtot = getVtot(mol, vion, vext, vhart, vxc);
        H.vtot = vtot;
        H.rho = rho;
    end
    Etotvec(iterdcm) = real(Etot);
    % update the total potential
    P = Y(:,nocc+1:end)*G(nocc+1:end,1:nocc);
    
    vtoterr = norm(vtot(:)-vtotin(:))/norm(vtotin(:));
    rhoerr = norm(rho(:)-rhoin(:))/norm(rhoin(:));
    trdcmerr(iterdcm) = max(vtoterr,rhoerr);
    fprintf('Rel Vtot Err    = %20.3e\n',vtoterr);
    
    % Convergence check
    fprintf('Total Energy    = %20.13e\n', Etot);
    [cvg,resfro] = reportconverge(H,X,iterdcm,maxdcmiter, ...
        trdcmerr(iterdcm),trdcmtol,verbose);
    if cvg
        info.converge = true;
        break;
    end
end

X.occ = zeros(1,ncols(X));
X.occ(1:nocc) = 1;
if force
    mol.xyzforce = getFtot(mol,H,X,rho);
end

info.Eigvals = ev;
info.Etotvec = Etotvec(1:iterdcm);
info.TRDCMerrvec = trdcmerr(1:iterdcm);

timetot = cputime - tstart;
fprintf('Etot            = %20.13e\n', Etot);
fprintf('Ekin            = %20.13e\n', Ekin);
fprintf('Eewald          = %20.13e\n', Eewald);
fprintf('Ealphat         = %20.13e\n', Ealphat);
fprintf('Ecor            = %20.13e\n', Ecor);
fprintf('Ecoul           = %20.13e\n', Ecoul);
fprintf('Exc             = %20.13e\n', Exc);
fprintf('--------------------------------------\n');
fprintf('Total time used = %20.3e\n', timetot);
fprintf('||HX-XD||_F     = %20.3e\n', resfro);

end
