function w = ncols(BX)
% BLOCHWAVEFUN/NCOLS Numbers of columns in the Bloch wave function class
%    w = NCOLS(BX) returns the numbers of columns in the Bloch wave
%    function.
%
%    See also BlochWavefun.

%  Copyright (c) 2016-2017 Yingzhou Li and Chao Yang,
%                          Stanford University and Lawrence Berkeley
%                          National Laboratory
%  This file is distributed under the terms of the MIT License.

w = cellfun(@(X)ncols(X),BX.wavefuncell);

end
