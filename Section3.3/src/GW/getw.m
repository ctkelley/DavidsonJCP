function W = getw(ksinfo, chi0)
%
%  W = getw(ksinfo, chi0);
%
%  compute the static screened Coulomb from irreducible polarizability chi0
%  and the Coulomb interaction in reciprocal space
%
ng = size(chi0,1);
vol  = ksinfo.vol;
epsilon = geteps(chi0, ksinfo.coulG);
ksinfo.coulG(1) = ksinfo.coulG0;
Dcoul   = spdiags(ksinfo.coulG,0,ng,ng);
%W = inveps*Dcoul/vol;
%inveps = inveps(2:end,2:end);
%Dcoul = Dcoul(2:end,2:end);
W = epsilon\Dcoul;
