function [X,eps,coulG,rsgrid,kdim,kpts,Nv,VA] = test_isdfk_gen()
% Generates the data needed for testing ISDF with k-points.
% OUTPUT
%   X      -- cell array of wavefunctions.  each corresponds to a different
%             k-point.  (output from test_isdfk_gen.m)
%   eps    -- matrix of size (Nbands)x(Nk) that holds the eigenvalues
%   coulG  -- vector that represents the diagonal of the Coulomb operators
%             in reciprocal space
%   rsgrid -- structure that contains the real space grid points in the 
%             x, y, z directions.  grid.x, grid.y, grid.z
%   kdim   -- dimensions of the k-point grid
%   kpts   -- the k-points
%   Nv     -- number of valence (occupied) orbitals
%   VA     -- the exchange kernel V_A in the BSE

    % make sure starting from correct directory
    try
        run('../../KSSOLV_startup.m');
    catch err
        error('Must start from .../crd-kssolv/src/ISDF');
    end

    % ---------------------------------------------------------------------
    % Define parameters
    % ---------------------------------------------------------------------
    Nbands = 14; % total number of bands (Nc + Nv)
    
    % There is an "autokpts" feature in sih4_scf_k.m that can be used to
    % include more kpts.
    
    % ---------------------------------------------------------------------
    % Construct crystal and do scf
    % ---------------------------------------------------------------------
    run('../../example/si2_scf_k.m');
    % This outputs:
    % [cry,H,X,info] = scf(cry);
    % where H is BlochHam
    %       X is BlochWavefun
    wf = X{1};
    Nk = H.nkpts;              % number of k points
    Nv = size(wf.occ,1);       % number of valence (occ) bands
    Nc = Nbands - Nv;          % number of conduction bands
    Ngrid = H.n1*H.n2*H.n3;    % number of spatial grid points
    Nfourier = size(wf.psi,1); % number of Fourier grid points
    clear X wf;
    
    % ---------------------------------------------------------------------
    % Get conduction/valence bands at all k-points
    % ---------------------------------------------------------------------
	X = cell(Nk,1);
    eps = zeros(Nbands,Nk);
    
    for ik=1:Nk
        disp(horzcat('Getting eigenvectors for k-point: ' ...
                    ,num2str(ik),'/',int2str(Nk)));
        [X{ik},eps(:,ik)] = diagbyeigs(cry,H{ik},Nbands,10^-10,300);
    end
    
    % ---------------------------------------------------------------------
    % Construct the bare Coulomb in reciprocal space
    % ---------------------------------------------------------------------
    fprintf('Computing Coulomb... ');
    grid = Ggrid(cry);
    gkk = grid.gkk; % get a vector |G|^2 within the Ecut limit
    idxnz = grid.idxnz; % get the position of gkk within a cube
    coulG = zeros(Nfourier,1);
    scell = cry.supercell;
    amin = min([norm(scell(:,1)) norm(scell(:,2)) norm(scell(:,3))])/2; % not quite right
    for j = 1:length(idxnz)
       if ( abs(gkk(j)) ~= 0 )
          coulG(j) = 4.0*pi/gkk(j); 
          coulG(j) = coulG(j)*(1-cos(sqrt(gkk(j))*amin));
       else
          %coulG(j) = 4.0*pi*amin^2/2;
          coulG(j) = 0.0;
       end
    end
    fprintf('Done!\n');
    
    % ---------------------------------------------------------------------
    % Determine dimensions of k-grid
    % ---------------------------------------------------------------------
    kdim = zeros(3,1);
    kpts = cry.kpts;
    k1k2k3 = size(kpts,1);
    fk = find(kpts(:,3) ~= kpts(1,3));
    if isempty(fk)
        kdim(3) = 1;
        k1k2 = k1k2k3;
    else
        k1k2 = fk(1) - 1;
    end
    fk = find(kpts(:,2) ~= kpts(1,2));
    if isempty(fk)
        kdim(2) = 1;
        kdim(1) = k1k2;
    else
        kdim(1) = fk(1) - 1;
    end
    kdim(2) = k1k2/kdim(1);
    kdim(3) = k1k2k3/k1k2;
    
    % ---------------------------------------------------------------------
    % Construct real space grid
    % ---------------------------------------------------------------------
    
    supercell = cry.supercell;
    if(nnz(supercell([2,3,4,6,7,8])) ~= 0)
        error('Only diagonal supercells supported.');
    end
    n1 = cry.n1;
    n2 = cry.n2;
    n3 = cry.n3;
    rsgrid.x = reshape(linspace(0, (n1-1)/n1*supercell(1,1), n1), [], 1, 1);
    rsgrid.y = reshape(linspace(0, (n2-1)/n2*supercell(2,2), n2), 1, [], 1);
    rsgrid.z = reshape(linspace(0, (n3-1)/n3*supercell(3,3), n3), 1, 1, []);
    
    % ---------------------------------------------------------------------
    % Construct W
    % ---------------------------------------------------------------------
    
%     cd ../GW;
%     ksinfo = gwsetup(cry);
%     omega = 0;
%     eta = 0;
%     [ng, nr] = size(ksinfo.F);
%     nclimit = 300;
%     chi0 = getchi0(ksinfo, omega, eta, nclimit);
%     epsilon = geteps(chi0, ksinfo.coulG);
%     inveps = inv(epsilon);
%     W = getw(ksinfo, chi0);
    
    % ---------------------------------------------------------------------
    % Construct exchange kernel: V_A(v,c,k,v',c',k')
    % ---------------------------------------------------------------------
    
    % Transform X to spatial basis
    starttrans = tic;
    Xrs = zeros(Ngrid, Nbands, Nk);
    for ik = 1:Nk
        wftemp = ifft3(X{ik});
        Xrs(:,:,ik) = wftemp.psi;
    end
    fprintf('Transform to spatial basis: %f sec\n', toc(starttrans));
    clear wftemp;
    
    % Construct orbital pairs: rho_{vck}(x) = psi_{vk}^*(x) psi_{ck}(x)
    startrho = tic;
    rho = cell(Nk,1);
    for ik = 1:Nk
        rho{ik} = Wavefun();
        rho{ik}.n1 = n1; rho{ik}.n2 = n2; rho{ik}.n3 = n3;
        rho{ik}.psi = zeros(Ngrid,Nv,Nc);
        for ic = 1:Nc
            rho{ik}.psi(:,:,ic) = ...
                conj(Xrs(:,1:Nv,ik)).*repmat(Xrs(:,Nv+ic,ik),1,Nv);
        end
        rho{ik}.psi = reshape(rho{ik}.psi, Ngrid, Nv*Nc);
    end
    fprintf('Construct rho: %f sec\n', toc(startrho));
    
    % Transform rho to Fourier basis
    starttrans = tic;
    rhof = zeros(Nfourier, Nv, Nc, Nk);
    for ik = 1:Nk
        tempwf = fft3(rho{ik});
        tempwf = tempwf(idxnz,:);
        rhof(:,:,:,ik) = reshape(tempwf.psi, Nfourier, Nv, Nc);
    end
    rhof = reshape(rhof, Nfourier, Nv*Nc*Nk);
    fprintf('Transform to Fourier basis: %f sec\n', toc(starttrans));
    clear rho tempwf;
    
    % Construct exchange kernel
    startExact = tic;
    fprintf('Calculating exchange kernel... ');
    VA = rhof'*spdiags(coulG,[0],length(coulG),length(coulG))*rhof;
    fprintf('%f seconds.\n\n',toc(startExact));
    clear rhof;
    
end
