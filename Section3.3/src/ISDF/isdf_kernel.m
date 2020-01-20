function [zeta_mu] = isdf_kernel(Phi, Psi, ind_mu, rk)

randn('state', 0);

[m1, n1] = size(Phi);
[m2, n2] = size(Psi);

if m1 ~= m2,
    error('Wrong inputs: row dimensions of Phi and Psi do not match!\n');
end
m = m1;
if length(ind_mu) > m,
    error('Wrong inputs: ind_mu is too long!\n');
end

%G = speye(n1*n2);
G = randn(n1*n2, ceil(2.0*rk));
%G = eye(n1*n2)
B = prod_states(Phi, Psi);
BG = B*G;
AG = BG(ind_mu, :);
zeta_mu = BG*pinv(AG);

return;
