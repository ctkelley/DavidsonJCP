function Ecor = getEcor(mol, rho, vtot, vion, vext)
%
% Usage: Ecor = getEcor(mol, rho, vtot, vion, vext)
%
% Purpose:
%    Computes an energy correction term
%
% Input:
%    mol  --- a Molecule object
%    rho  --- charge density (3D array)
%    vtot --- total local potential (3D array)
%    vion --- local ionic potential (3D array)
%    vext --- external potential (3D array)
%
% Ouptut:
%    Ecor --- correction energy
%
hhh = mol.vol/mol.n1/mol.n2/mol.n3;
Ecor      = sumel( (vion+vext-vtot).*rho );
Ecor      = real(Ecor)*hhh;
