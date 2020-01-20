function [mol,H,X,info] = dcm(mol,options)
% DCM Direct Constrained Minimization algorithm.
%    [mol,H,X,info] = DCM(mol,options) adopts Direct Constrained
%    Minimization algorithm to find the ground state minimum total energy
%    and the corresponding wave functions. mol is a Molecule object and
%    options is the options for running the DCM. Please read setksopt for
%    detailed information about options. DCM returns the molecule mol, the
%    Hamiltonian H, the wave functions X, and the information about
%    iteractions. If a crystal object is input instead of molecule, the
%    results will be generated by dcm4c function.
%
%   See also scf, trdcm, dcm4m, dcm4c.

%  Copyright (c) 2016-2017 Yingzhou Li and Chao Yang,
%                          Stanford University and Lawrence Berkeley
%                          National Laboratory
%  This file is distributed under the terms of the MIT License.

if (nargin < 2)
    options = setksopt();
end

if mol.temperature > 0
    [mol,H,X,info] = dcm4m(mol,options);
    return;
elseif isa(mol,'Crystal')
    [mol,H,X,info] = dcm4c(mol,options);
    return;
end

% Set timer
tstart  = cputime;

% Initialize input variables
verbose    = ~strcmpi(options.verbose,'off');
maxdcmiter = options.maxdcmiter;
maxinerscf = options.maxinerscf;
dcmtol     = options.dcmtol;
X          = options.X0;
rho        = options.rho0;

nspin      = mol.nspin;

% Initialize Hamiltonian, Wavefun, and Preconditioners
[mol,H,X,Hprec,nocc] = iterinit(mol,rho,X);

% calculate Ewald and Ealphat
Eewald     = getEewald(mol);
Ealphat    = getEalphat(mol);

vion       = H.vion;
vext       = H.vext;

% Initialize output variables
Etotvec    = zeros(maxdcmiter,1);
dcmerr     = zeros(maxdcmiter,1);

fprintf('Beging DCM calculation for %s...\n',mol.name);

[X,~] = updateX(mol, H, X, Hprec, options);
rho = getcharge(mol,X);
[vhart,vxc,uxc2,rho]=getVhxc(mol,rho);
vtot = getVtot(mol, vion, vext, vhart, vxc);
H.vtot = vtot;

for iterdcm = 1:maxdcmiter
    fprintf('DCM iter %2d:\n', iterdcm);
    
    vtotin = vtot;
    
    % calculate the residual
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
    
    % solve the projected problem
    B = real(Y'*Y); B = (B+B')/2;
    for iterscf = 1:maxinerscf
        A = real(Y'*(H*Y));
        A = (A+A')/2;
        
        [G,D]=eig(A,B,'chol');
        G = G(:,1:nocc);
        ev = diag(D);
        
        X = Y*G(:,1:nocc);
        rho = getcharge(mol,X);
        
        % Kinetic energy and some additional energy terms
        Ekin = (2/nspin)*sum(ev(1:nocc));
        % ionic and external potential energy was included in Ekin
        % along with incorrect Ecoul and Exc. Need to correct them
        % later;
        Ecor = getEcor(mol, rho, vtot, vion, vext);
        
        % Compute Hartree and exchange correlation energy and potential
        % using the new charge density; update the total potential
        [vhart,vxc,uxc2,rho]=getVhxc(mol,rho);
        
        % Update total potential
        vtot = getVtot(mol, vion, vext, vhart, vxc);
        % TODO: Should update vnp instead of vtot, at the same time, vtot
        % should be consist with the vnp
        H.vtot = vtot;
        H.rho = rho;
        
    end
    
    % Calculate the potential energy based on the new potential
    Ecoul = getEcoul(mol,abs(rho),vhart);
    Exc   = getExc(mol,abs(rho),uxc2);
    Etot = Eewald + Ealphat + Ekin + Ecor + Ecoul + Exc;
    Etotvec(iterdcm) = Etot;
    
    % update the total potential
    P = Y(:,nocc+1:end)*G(nocc+1:end,:);
    
    vtoterr = norm(vtot(:)-vtotin(:))/norm(vtotin(:));
    dcmerr(iterdcm) = vtoterr;
    fprintf('Rel Vtot Err    = %20.3e\n',vtoterr);

    % Convergence check
    fprintf('Total Energy    = %20.13e\n', Etot);
    [cvg,resfro] = reportconverge(H,X,iterdcm,maxdcmiter, ...
        dcmerr(iterdcm),dcmtol,verbose);
    if cvg
        info.converge = true;
        break;
    end
end

X.occ = zeros(1,ncols(X));
X.occ(1:nocc) = 1;
mol.xyzforce = getFtot(mol,H,X,rho);

info.Eigvals = ev;
info.Etotvec = Etotvec(1:iterdcm);
info.DCMerrvec = dcmerr(1:iterdcm);

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
