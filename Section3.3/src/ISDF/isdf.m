function [zeta_mu, ind_mu] = isdf(Phi, Psi, rk)
% The function isdf() performs the
% Interpolative Separable Density Fitting
% decomposition:
%   phi_i(r)*psi_j(r) = sum_{mu} zeta_mu(r)*phi_i(r_mu)*psi_j(r_mu).
%
% Inputs:
%   Phi: A set of wavefunctions (size: M*N1).
%   Psi: Another set of wavefunctions (size: M*N2).
%   rk: rank of product states [ Phi_i*Psi_j ],
%       1 <= rk <= min(M, N1*N2).
% Outputs:
%   zeta_mu: interpolation vectors (basis) for density fitting
%   ind_mu: indices of interpolation points

% Find the interpolation points.
ind_mu = samp(Phi, Psi, rk);

% Construct the interpolation basis.
zeta_mu = isdf_kernel(Phi, Psi, ind_mu, rk);

return;
