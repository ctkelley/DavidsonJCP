function B = repcol(A,ncols)
% REPCOL Replicate each columns of a matrix.
%     B = REPCOL(A,ncols) returns a matrix B with the same number of rows
%     as A and number of columns as sum(ncols). The ith column of A is
%     replicated ncols(i) times contigously.
%
%   Example:
%     >> A = [ 1 2
%              3 4
%              5 6 ];
%     >> B = repcol(A,[2,3])
%     
%     B =
%           1 1 2 2 2
%           3 3 4 4 4
%           5 5 6 6 6
%
%   Remark: This function currently is not well optimized. Please let the
%   authors know if there is any better implementation.
%  
%     See also repmat.

%  Copyright (c) 2015-2016 Yingzhou Li and Chao Yang,
%                          Stanford University and Lawrence Berkeley
%                          National Laboratory
%  This file is distributed under the terms of the MIT License.

Acell = mat2cell(A,size(A,1),ones(1,size(A,2)));
Acell = cellfun(@(A,b)repmat(A,1,b),Acell,num2cell(ncols), ...
    'UniformOutput',0);
B = cell2mat(Acell);

end