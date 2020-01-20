function ind_mu = samp(Phi, Psi, rk)

rand('state', 0);

[m1, n1] = size(Phi);
[m2, n2] = size(Psi);

if m1 ~= m2,
    error('Wrong inputs: row dimensions of Phi and Psi do not match\n');
end
m = m1;

% Compute product states
B = prod_states(Phi, Psi);

% Do dimension reduction via random Guassian projection
G = randn(n1*n2, ceil(2*rk));
BG = B*G;
clear B G;

% Do QRCP
[~, ~, e] = qr(BG',0);

% Indices of important columns
ind_mu = e(1:rk);

return
