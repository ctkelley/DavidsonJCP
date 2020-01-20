function Fewald = getFewald(mol)

nel = mol.nel;
alist = mol.alist;
atoms = mol.atoms;
xyzlist = mol.xyzlist;
na = sum(mol.natoms);
ecut2    = mol.ecut2;
e2 = e2Def();

grid2 = Ggrid(mol,ecut2);
gkk   = grid2.gkk;
gkx   = grid2.gkx;
gky   = grid2.gky;
gkz   = grid2.gkz;

alpha = 1.1;
upperbound = 1;
while upperbound >1e-6
    alpha = alpha-0.1;
    upperbound = e2*nel^2*sqrt(2*alpha/2/pi)*erfc(sqrt(ecut2/2/alpha));
end

if abs(gkk(1)) <= 1e-8
    idxg = 2:length(gkk);
else
    idxg = 1:length(gkk);
end

strf = exp( -1i*[gkx(idxg) gky(idxg) gkz(idxg)]*xyzlist' );
aux = sum( conj(strf)*diag(arrayfun(@(x)atoms(x).venum,alist) ), 2 ) ...
    .*exp(-gkk(idxg)/alpha/4)./gkk(idxg);

Fewald = zeros(na,3);

for it = 1:na
    arg = [gkx(idxg) gky(idxg) gkz(idxg)]*xyzlist(it,:)';
    sumnb = cos(arg).*imag(aux) - sin(arg).*real(aux);
    fact = -atoms(alist(it)).venum*e2*4*pi/mol.vol;
    Fewald(it,1) = fact*(gkx(idxg)'*sumnb);
    Fewald(it,2) = fact*(gky(idxg)'*sumnb);
    Fewald(it,3) = fact*(gkz(idxg)'*sumnb);
end

if abs(gkk(1)) > 1e-8
    return;
end

rmax = 5/sqrt(alpha);

for iti = 1:na
    for itj = 1:na
        if iti ~= itj
            dxyz = xyzlist(iti,:)-xyzlist(itj,:);
            [r,r2] = rgen(dxyz,rmax);
            if size(r,1) == 0
                continue;
            end
            rr = sqrt(r2);
            factor = atoms(alist(iti)).venum*atoms(alist(itj)).venum ...
                *e2./rr.^2.*( erfc(sqrt(alpha)*rr)./rr ...
                + sqrt(8*alpha/2/pi)*exp(-alpha*rr.^2) );
            Fewald(iti,1) = Fewald(iti,1) - r(:,1)'*factor;
            Fewald(iti,2) = Fewald(iti,2) - r(:,2)'*factor;
            Fewald(iti,3) = Fewald(iti,3) - r(:,3)'*factor;
        end
    end
end

    function [r,r2] = rgen(dxyz,rmax)
        at = mol.supercell;
        bg = inv(mol.supercell);
        n1 = 2*floor(norm(bg(:,1))*rmax)+2;
        n2 = 2*floor(norm(bg(:,2))*rmax)+2;
        n3 = 2*floor(norm(bg(:,3))*rmax)+2;
        
        ds = dxyz*bg;
        ds = ds-round(ds);
        dxyz = ds*at;
        
        [I,J,K] = ndgrid((0:n1-1)-((0:n1-1) > n1/2)*n1, ...
                    (0:n2-1)-((0:n2-1) > n2/2)*n2, ...
                    (0:n3-1)-((0:n3-1) > n3/2)*n3);
                
        r = [I(:),J(:),K(:)]*at'-repmat(dxyz,n1*n2*n3,1);
        r2 = sum(r.^2,2);
        idx = r2 <= rmax^2 & r2 > 1e-10;
        r = r(idx,:);
        r2 = r2(idx);
        
    end

end
