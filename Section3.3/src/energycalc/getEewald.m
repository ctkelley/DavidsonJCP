function Eewald = getEewald(mol)

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

alpha = 2.9;
upperbound = 1;
while upperbound >1e-7
    alpha = alpha-0.1;
    upperbound = 2*nel^2*sqrt(alpha/pi)*erfc(sqrt(ecut2*2*meDef()/4/alpha));
end

% sum in G-space
if abs(gkk(1)) <= 1e-8
    idxg = 2:length(gkk);
    Eewaldg = -nel^2/alpha/4;
else
    idxg = 1:length(gkk);
    Eewaldg = 0;
end

strf = exp( -1i*[gkx(idxg) gky(idxg) gkz(idxg)]*xyzlist' );
aux = abs( sum( conj(strf) ...
    *diag(arrayfun(@(x)mol.ppvar.venums(x),alist) ), 2 )).^2 ...
    .*exp(-gkk(idxg)/alpha/4)./gkk(idxg);

Eewaldg = (Eewaldg + sum(aux))*2*pi/mol.vol;

if abs(gkk(1)) <= 1e-8
    Eewaldg = Eewaldg ...
        - sqrt(1/pi*alpha)*sum(arrayfun(@(x)mol.ppvar.venums(x),alist).^2);
end

% sum in R-space
Eewaldr = 0;
rmax = 4/sqrt(alpha);

for iti = 1:na
    for itj = 1:na
            dxyz = xyzlist(iti,:)-xyzlist(itj,:);
            [r,r2] = rgen(dxyz,rmax);
            if size(r,1) == 0
                continue;
            end
            rr = sqrt(r2);
            Eewaldr = Eewaldr + ...
                atoms(alist(iti)).venum*atoms(alist(itj)).venum ...
                *sum(erfc(sqrt(alpha)*rr)./rr);
    end
end

% TODO: figure out why divide by 2 or multiply Eewaldg by 2
Eewald = e2*(Eewaldg+Eewaldr/2);

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
