function Eion = getEion(mol,rho,vion)
%
% Usage: Eion = getEion(rho, vion);
%
% Purpose:
%    Compute the ionic potential energy
%
% Input:
%    mol  --- a Molecule object
%    rho  --- input charge density (3D array)
%    vion --- ionic potential (3D array)
%
% Output:
%   Eion  --- Ionic potential energy (scalar)
%
hhh = mol.vol/mol.n1/mol.n2/mol.n3;
Eion      = sumel(rho.*vion);
Eion      = real(Eion)*hhh;
