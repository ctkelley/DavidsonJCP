function Fnl = getFnl(mol,H,X)

na = sum(mol.natoms);
ecut = mol.ecut;

grid2 = Ggrid(mol,ecut);
gkx   = grid2.gkx;
gky   = grid2.gky;
gkz   = grid2.gkz;

[llsiz,lloffset] = getvnlsize(mol);

if isa(mol,'Crystal')
    nkpts = mol.nkpts;
    Fnl = zeros(na,3);
    for ik = 1:nkpts
        Fnl = Fnl + mol.wks(ik) ...
            *getFnlsingle(H.vnlmatcell{ik},H.vnlsigncell{ik},X{ik});
    end
else
    Fnl = getFnlsingle(H.vnlmat,H.vnlsign,X);
end

    function Fnlsingle = getFnlsingle(bec,sig,X)
        
        X = X*diag(X.occ);
        
        becX = bec'*X;
        
        dbecX = bec'*(repmat(1i*gkx,1,ncols(X)).*X);
        dbecXbecX = real(becX.*conj(dbecX));
        Fx = -4*sum(repmat(sig,1,size(becX,2)).*dbecXbecX,2);
        
        dbecX = bec'*(repmat(1i*gky,1,ncols(X)).*X);
        dbecXbecX = real(becX.*conj(dbecX));
        Fy = -4*sum(repmat(sig,1,size(becX,2)).*dbecXbecX,2);
        
        dbecX = bec'*(repmat(1i*gkz,1,ncols(X)).*X);
        dbecXbecX = real(becX.*conj(dbecX));
        Fz = -4*sum(repmat(sig,1,size(becX,2)).*dbecXbecX,2);
        
        Fnlsingle = zeros(na,3);
        
        for ia = 1:na
            Fnlsingle(ia,1) = sum(Fx(lloffset(ia)+(1:llsiz(ia))));
            Fnlsingle(ia,2) = sum(Fy(lloffset(ia)+(1:llsiz(ia))));
            Fnlsingle(ia,3) = sum(Fz(lloffset(ia)+(1:llsiz(ia))));
        end
    end

end
