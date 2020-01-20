function test_isdfk(X, coulG, Nv, rsgrid, kdim, kpts, VA, ...
                    doISDFon, isdfkcofac, Nu)
% This function performs ISDF and then checks it to see if it is accurate.
% INPUT
%   X      -- cell array of wavefunctions.  each corresponds to a different
%             k-point.  (output from test_isdfk_gen.m)
%   coulG  -- vector that represents the diagonal of the Coulomb operators
%             in reciprocal space.  (output from test_isdfk_gen.m)
%   Nv     -- number of valence orbitals
%   rsgrid -- structure that contains the real space grid points in the 
%             x, y, z directions.  grid.x, grid.y, grid.z
%   kdim   -- 3-vector that gives the number of k-points in each dimension
%             of the (fine) k-grid.
%   kpts   -- the k-points
%   VA     -- exhange kernel in BSE
%   doISDFon   -- 'psi' to do ISDF on psi.  Otherwise it will perform ISDF
%                 on u.
%   isdfkcofac -- the factor by which to reduce the k-point in *each*
%                 dimension when doing the ISDF.  isdfkcofac = 1 does not
%                 reduce the amount of k-points.  each dimension of 
%                 k-points must be a factor of isdfkcofac.  Example:
%                 kdim: [4,4,4], and isdfkcofac = 2.  then ISDf will use 
%                 the k-points for a 2x2x2 k-point grid.
%   Nu     -- number of auxiliary basis functions
    
    Nk = size(X,1);
    [Nfourier, Nbands] = size(X{1}.psi);
    n1 = X{1}.n1;
    n2 = X{1}.n2;
    n3 = X{1}.n3;
    idxnz = X{1}.idxnz;
    Ngrid = n1*n2*n3;
    Nc = Nbands - Nv;
    
    % Check if isdfkcofac is valid
    if sum(mod(kdim,isdfkcofac) ~= 0) ~= 0 || isdfkcofac < 1
        error('isdfkfac is invalid');
    end
    
    % Determine coarse k-point grid
    [x,y,z] = ndgrid(1:kdim(1), 1:kdim(2), 1:kdim(3));
    finegrid = [x(:), y(:), z(:)];
    inc = isdfkcofac;
    [x,y,z] = ndgrid(1:inc:kdim(1), 1:inc:kdim(2), 1:inc:kdim(3));
    coarsegrid = [x(:), y(:), z(:)];
    idxco = zeros(size(coarsegrid,1),1);
    for i = 1:size(coarsegrid,1)
        for j = 1:size(finegrid,1)
            if sum(coarsegrid(i,:) == finegrid(j,:)) == 3
                idxco(i) = j;
                break;
            end
            if j == size(finegrid,1)
                error('Did not find k-point'); % should never reach here
            end
        end
    end
    Nkco = length(idxco);
    clear i j x y z inc findgrid coarsegrid;
    
    % ---------------------------------------------------------------------
    % Construct exchange kernel via ISDF
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

    % Convert u to psi
    %   - Comment this part out if you want to use u.
    Xrs = reshape(Xrs, n1, n2, n3, Nbands, Nk);
    if strcmp(doISDFon,'psi')
        for ik = 1:Nk
            for ib = 1:Nbands
                Xrs(:, :, :, ib, ik) = exp(1i * ( ...
                    repmat(rsgrid.x, 1, n2, n3) * kpts(ik,1) ...
                    + repmat(rsgrid.y, n1, 1, n3) * kpts(ik,2) ...
                    + repmat(rsgrid.z, n1, n2, 1) * kpts(ik,3))) ...
                    .* Xrs(:, :, :, ib, ik);
            end
        end
    end
    Xrs = reshape(Xrs, Ngrid, Nbands, Nk);
    % END convert u to psi.
    
    Xco = reshape(Xrs(:,:,idxco), Ngrid, Nbands*Nkco); % on coarse grid
    Xrs = reshape(Xrs, Ngrid, Nbands*Nk);
    
    % Perform ISDF on *all* pairs
    [P, idx_mu] = isdf_fast(conj(Xco), Xco, Nu);
    clear Xco;
    Xmu = Xrs(idx_mu,:); % coeffs for ISDF
    Xmu = reshape(Xmu, Nu, Nbands, Nk);
    
    % Test ISDF
    %  - WARNING: This takes *forever* for even medium size systems.
%     Xrs = reshape(Xrs, Ngrid, Nbands, Nk);
%     rho = zeros(Ngrid, Nbands^2);
%     err2 = 0;
%     rhonorm2 = 0;
%     for ik1 = 1:Nk
%         for ik2 = 1:Nk
%             k = 1;
%             for i = 1:Nbands
%                 for j = 1:Nbands
%                     rho(:, k) = conj(Xrs(:, i, ik1)).*Xrs(:, j, ik2);
%                     k = k + 1;
%                 end
%             end
%             err2 = err2 + norm(P*rho(idx_mu, :) - rho, 'fro')^2;
%             rhonorm2 = rhonorm2 + norm(rho, 'fro')^2;
%         end
%     end
%     fprintf('Rel error in ISDF: %e\n\n', sqrt(err2/rhonorm2));
    clear Xrs rho;
    
    % Transform P to Fourier basis
    Pwf = Wavefun();
    Pwf.n1 = n1; Pwf.n2 = n2; Pwf.n3 = n3;
    Pwf.psi = P;
    Pwf = fft3(Pwf);
    Pwf = Pwf(idxnz,:);
    P = Pwf.psi;
    clear Pwf;
    
    % Construct <P_mu | v | P_nu>
    startV = tic;
    V = P'*spdiags(coulG, [0], length(coulG), length(coulG))*P;
    fprintf('Construct V: %f sec\n', toc(startV));
    clear P;
    
    % Construct C.  C(mu,v,c,k) = psi_{vk}^*(x_mu) psi_{ck}(x_mu)
    startC = tic;
    C = zeros(Nu, Nv, Nc, Nk);
    for ik = 1:Nk
        for ic = 1:Nc
            C(:,:,ic,ik) = ...
                conj(Xmu(:,1:Nv,ik)).*repmat(Xmu(:,Nv+ic,ik),1,Nv);
        end
    end
    fprintf('Construct C: %f sec\n', toc(startC));
    
    % Construct exchange kernel
    %   VA_isdf is of size (Nv*Nc*Nk) x (Nv*Nc*Nk) and is in the proper 
    %   order to be reshaped to Nv x Nc x Nk x Nc x Nv x Nk.
    startkernel = tic;
    C = reshape(C, Nu, Nv*Nc*Nk);
    VA_isdf = C'*V*C;
    fprintf('Construct ISDF kernel: %f sec\n', toc(startkernel));
    
    % Compare kernels
    VA = reshape(VA,Nv^2*Nc^2*Nk^2,1);
    VA_isdf = reshape(VA_isdf,Nv^2*Nc^2*Nk^2,1);
    fprintf('Relative error in V_A: %e\n\n', norm(VA-VA_isdf)/norm(VA));
    
end
