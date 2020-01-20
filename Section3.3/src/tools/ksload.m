function [mol,H,X,info] = ksload(savefile)
% KSLOAD load KSSOLV related data into workspace.
%
%   [mol,H,X,info] = KSLOAD(savefile) loads molecule, hamiltonian,
%   wavefunction and scf iteration information from a file named as
%   savefile.
%
%   See also kssave.

%  Copyright (c) 2016-2018 KSSOLV Team
%  This file is distributed under the terms of the MIT License.

ksdat = load(savefile);
mol = ksdat.mol;
H = ksdat.H;
X = ksdat.X;
info = ksdat.info;

end