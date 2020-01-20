function [vnlmat,vnlsign] = getvnl(mol,kpts)
% GETVNL returns the non-local pseudo potential of the molecule.
%    [vnlmat,vnlsign] = GETVNL(mol) calculate the non-local pseudo
%    potential of the molecule based on the non-local pseudo potential of
%    each atom in the molecule in G-space.
%
%    See also vnl2g, vloc2g, getvloc.

%  Copyright (c) 2015-2016 Yingzhou Li and Chao Yang,
%                          Stanford University and Lawrence Berkeley
%                          National Laboratory
%  This file is distributed under the terms of the MIT License.

atoms = mol.atoms;
alist = mol.alist;
ntypes = length(atoms);
na = length(alist);
vol = mol.vol;
dq = 0.01;

if nargin == 1
    kpts = [ 0 0 0 ];
end

ppvar = mol.ppvar;

xyzlist  = mol.xyzlist;

grid = Ggrid(mol);
ng    = grid.ng;
gkx   = grid.gkx+kpts(1);
gky   = grid.gky+kpts(2);
gkz   = grid.gkz+kpts(3);
kkxyz = [gkx gky gkz];
gkk   = sum(kkxyz.*kkxyz,2);
gk    = sqrt(gkk);
%
phmat  = [gkx gky gkz]*xyzlist';
sk = exp(1i.*phmat);
vnlcell = cell(1,na);
vnlD = cell(1,na);
maxl = max(cell2mat(ppvar.lll));
ylmrtable = ylmr(gkx,gky,gkz,maxl);

% atypes should list all types ordered by atomic numbers
for itype = 1:ntypes
    nbeta = size(ppvar.vnltab{itype},2);
    if nbeta == 0
        index  = find(alist == itype);
        for ita = index'
            vnlcell{ita} = zeros(ng,0);
            vnlD{ita} = [];
        end
        continue;
    end
    totall = 0;
    for it = 1:length(ppvar.lll{itype})
        ll = ppvar.lll{itype}(it);
        totall = totall + 2*ll+1;
    end
    vkb1 = zeros(ng,totall);
    tab = ppvar.vnltab{itype};
    px = repmat(gk/dq-round(gk/dq),1,nbeta);
    ux = 1-px;
    vx = 2-px;
    wx = 3-px;
    it0 = round(gk/dq)+1;
    it1 = it0+1;
    it2 = it0+2;
    it3 = it0+3;
    vq = tab(it0,:).*ux.*vx.*wx/6 ...
        + tab(it1,:).*px.*vx.*wx/2 ...
        - tab(it2,:).*px.*ux.*wx/2 ...
        + tab(it3,:).*px.*ux.*vx/6;
    offll = 0;
    for it = 1:length(ppvar.lll{itype})
        ll = ppvar.lll{itype}(it);
        idx = offll + (1:2*ll+1);
        tidx = ll^2 + (1:2*ll+1);
        vkb1(:,idx) = (-1i)^ll*ylmrtable(:,tidx) ...
            .*repmat(vq(:,it),1,2*ll+1);
        offll = offll + 2*ll + 1;
    end
    
    index  = find(alist == itype);
    for ita = index'
        vnlcell{ita} = repmat(conj(sk(:,ita)),1,totall).*vkb1;
    end
    
    D = ppvar.vnlD{itype};
    for ita = index'
        % The prefactor came from the Bessel transform in vkb and we
        % combined them here.
        vnlD{ita} = D*(4*pi)^2/vol;
    end
end

len = 0;
for it = 1:na
    len = len + size(vnlcell{it},2);
end
vnlmat  = zeros(length(gk),len);
vnlsign = zeros(len,1);
offset = 0;
for it = 1:na
    idx = offset + (1:size(vnlcell{it},2));
    [V,D] = eig(vnlD{it});
    vnlmat(:,idx) = vnlcell{it}*V;
    vnlsign(idx) = diag(D);
    offset = offset + size(vnlcell{it},2);
end

end
