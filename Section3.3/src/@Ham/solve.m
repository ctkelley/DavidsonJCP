function X = solve(H,B, shift)
% function to solve (H - shift*I)x = b
% where H is a hamiltonian

if nvargin == 2
    shift = zeros(1,ncol(B));
end


n1 = H.n1;
n2 = H.n2;
n3 = H.n3;
idxnz = H.idxnz;

xArray = zeros(length(idxnz), ncol(B));

for ii=1:ncol(B)
    % solving the problems using GMRES
    xArray(:,ii) = gmres( @(p) shift(ii)*p-mult(H,p), B.psi, [], 1e-6, 300 );
    
end
% Building the resultign Wavefield
X = Wavefun(xArray,n1,n2,n3,idxnz);

end