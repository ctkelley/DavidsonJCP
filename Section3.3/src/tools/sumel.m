function s = sumel(A)
% SUMEL sums all the elements.
%
%   s = SUMEL(A) returns the total sum of all elementals in A.
%
%   See also sum.

%  Copyright (c) 2016-2017 Yingzhou Li and Chao Yang,
%                          Stanford University and Lawrence Berkeley
%                          National Laboratory
%  This file is distributed under the terms of the MIT License.

if iscell(A)
    s = sumel(cellfun(@sumel,A));
else
    s = sum(A(:));
end

end