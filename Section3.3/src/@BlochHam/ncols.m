function w = ncols(BH)
% BLOCHHAM/NCOLS Numbers of columns in the Bloch Hamiltonian class
%    w = NCOLS(BH) returns the numbers of columns in the non-local pseudo
%    potential of the Bloch Hamiltonian.
%
%    See also BlochHam.

%  Copyright (c) 2016-2017 Yingzhou Li and Chao Yang,
%                          Stanford University and Lawrence Berkeley
%                          National Laboratory
%  This file is distributed under the terms of the MIT License.

w = cellfun(@(X)size(X,2),BH.vnlmatcell);

end
