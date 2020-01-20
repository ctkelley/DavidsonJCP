function m = minel(A)
% MINEL minimum of all the elements.
%
%   m = MINEL(A) returns the minimum of all elementals in A.
%
%   See also sumel, maxel.

%  Copyright (c) 2016-2017 Yingzhou Li and Chao Yang,
%                          Stanford University and Lawrence Berkeley
%                          National Laboratory
%  This file is distributed under the terms of the MIT License.

if iscell(A)
    m = minel(cellfun(@minel,A));
else
    m = min(A(:));
end

end