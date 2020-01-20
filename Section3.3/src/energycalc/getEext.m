function Eext = getEext(mol,rho,vext)
%
% Usage:
%    Eext = getEext(rho, vext);
%
% Purpose:
%    Compute the potential energy induced by an external potential vext.
%
% Input:
%  mol  --- Molecule object 
%  rho  --- input charge density (3D array)
%  vext --- external potential (3D array)
%
% Output:
%  Eext --- Potential energy induced by vext (scalar)
%

hhh = mol.vol/mol.n1/mol.n2/mol.n3;
Eext      = sumel(rho.*vext);
Eext      = real(Eext)*hhh;
