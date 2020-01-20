function P = prod_states(Phi, Psi)

[m1, n1] = size(Phi);
[m2, n2] = size(Psi);

if m1 ~= m2,
    error('Wrong inputs: row dimensions of Phi and Psi do not match!\n');
end
m = m1;

P = zeros(m, n1*n2);
for i2 = 1:n2,
    for i1 = 1:n1,
        P(:, i1 + (i2-1)*n1) = Phi(:, i1).*Psi(:, i2);
    end
end

return
