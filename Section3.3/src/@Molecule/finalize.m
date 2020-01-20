function mol = finalize(mol)
% MOLECULE/FINALIZE Finalize function for molecule class
%    mol = FINALIZE(mol) returns a molecule class of with finalized fields.
%
%    See also Molecule.

%  Copyright (c) 2016-2017 Yingzhou Li and Chao Yang,
%                          Stanford University and Lawrence Berkeley
%                          National Laboratory
%  This file is distributed under the terms of the MIT License.

if isempty(mol.ecut2)
    mol.ecut2 = 4*mol.ecut;
end

if isempty(mol.n1) || isempty(mol.n2) || isempty(mol.n3)
    C = mol.supercell;
    fackpt = sqrt(mol.ecut2*(2*meDef()))/pi;
    mol.n1 = ceil(fackpt*norm(C(1,:)));
    mol.n2 = ceil(fackpt*norm(C(2,:)));
    mol.n3 = ceil(fackpt*norm(C(3,:)));
end

if size(mol.atoms,1) < size(mol.atoms,2)
    mol.atoms = mol.atoms(:);
end

if size(mol.alist,1) < size(mol.alist,2)
    mol.alist = mol.alist(:);
end

if isempty(mol.vext)
    mol.vext = zeros(mol.n1,mol.n2,mol.n3);
end

if isempty(mol.nspin)
    mol.nspin = 1;
end

if isempty(mol.temperature)
    mol.temperature = 0.0;
end

if ~isempty(mol.atoms)
    mol.natoms = sum(repmat(mol.alist,1,length(mol.atoms)) == ...
        repmat(unique(mol.alist)',length(mol.alist),1),1);
end

mol.vol = abs(det(mol.supercell));

if isempty(mol.nel) && ~isempty(mol.atoms)
    mol.nel = sum(mol.natoms.*[mol.atoms.venum]);
end

end
