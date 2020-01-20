function [zeta_mu, ind_mu] = isdf_fast(Phi, Psi, rk)
% The function isdf_fast() performs the
% Interpolative Separable Density Fitting
% decomposition:
%   phi_i(r)*psi_j(r) = sum_{mu} zeta_mu(r)*phi_i(r_mu)*psi_j(r_mu).
% This function performs both steps of ISDF at the same time and is
% therefore capable of saving significant computational time in the second
% step by using the QR factorization from the first step to solve the least
% squares problem.  This function also uses the product structure to reduce
% the computation time of the random projection.  In this function, over
% 95% of the computation time is spent on the QRCP.
%
% Notation:
%   M: number of grid points
%   N1: number of Phi wavefunctions
%   N2: number of Psi wavefunctions
%
% Inputs:
%   Phi: A set of wavefunctions (size: M*N1).
%   Psi: Another set of wavefunctions (size: M*N2).
%   rk: rank of product states [ Phi_i*Psi_j ],
%       1 <= rk <= min(M, N1*N2).
% Outputs:
%   zeta_mu: interpolation vectors (basis) for density fitting
%   ind_mu: indices of interpolation points

    rng('default');

    [m1, n1] = size(Phi);
    [m2, n2] = size(Psi);

    if m1 ~= m2,
        error('Wrong inputs: row dimensions of Phi and Psi do not match\n');
    end
    startISDF = tic;

    % ---------------------------------------------------------------------
    % Find the interpolation points.
    % ---------------------------------------------------------------------

    % Do dimension reduction on Phi and Psi via random Guassian projection
    rsamp = 1.5*rk;
    r1 = min( ceil(sqrt((n1/n2)*rsamp)), n1 );
    r2 = min( ceil(sqrt((n2/n1)*rsamp)), n2 );
    G1 = randn(n1,r1);
    G2 = randn(n2,r2);
    PhiG = Phi*G1;
    PsiG = Psi*G2;
    
    % Compute product of the projected states
    BG = prod_states(PhiG, PsiG);

    % Do QRCP
    fprintf('Performing QRCP of size %d x %d...\n', size(BG,2), size(BG,1));
    startQRCP = tic;
    [~, R, e] = qr(BG',0);
    timeQRCP = toc(startQRCP);
    clear BG;

    % Indices of important columns
    ind_mu = e(1:rk);

    % ---------------------------------------------------------------------
    % Construct the interpolation basis.
    % ---------------------------------------------------------------------
    zeta_mu(e,:) = [eye(rk,rk), R(1:rk,1:rk)\R(1:rk,rk+1:end)]';
    timeISDF = toc(startISDF);
    
    % Print timing
    fprintf('QRCP: %f sec    ISDF: %f sec\n', ...
        timeQRCP, timeISDF);
    fprintf('%3.1f%% of ISDF time was spent on QRCP.\n', ...
        timeQRCP/timeISDF*100);
    
end