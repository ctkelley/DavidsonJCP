function [vhart,vxc,uxc,rho] = getVhxc(mol,rho)

% computing the potential from the Hartree approximation
vhart = getVhart(mol,rho);

% computing the exchange correlation, for rho and the molecule
[vxc,uxc,rho] = getVxc(mol,rho);

end