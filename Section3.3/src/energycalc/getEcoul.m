function Ecoul = getEcoul(mol,rho,vhart)
%
% Usage: Ecoul = getEcoul(mol, rho, vhart);
%
% Purpose:
%    Compute the Coulomb potential energy induced by vhart
%
% Input:
%  mol    --- a Molecule object
%  rho    --- input charge density (3D array)
%  vhart  --- Hartee (Coulomb) potential (3D array)
%
% Output:
%  Ecoul  --- Coulomb potential energy induced by vhart.
%
hhh = mol.vol/mol.n1/mol.n2/mol.n3;
Ecoul      = sumel(rho.*vhart)/2;
Ecoul      = real(Ecoul)*hhh;
