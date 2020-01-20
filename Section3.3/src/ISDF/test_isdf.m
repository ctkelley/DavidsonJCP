function err = test_isdf(psi, phi)

ng = size(psi, 1);
nv = size(psi, 2);
nc = size(phi, 2);
rk = min(ng, nv*nc);
[zeta_mu, ind_mu] = isdf(psi, conj(phi), rk);
V = [];
k = 1;
for i = 1:nv,
    for j = 1:nc,
        V(:, k) = psi(:, i).*conj(phi(:, j));
        k = k + 1;
    end
end
err = norm(zeta_mu*V(ind_mu, :) - V, 'fro')/norm(V, 'fro');

return
