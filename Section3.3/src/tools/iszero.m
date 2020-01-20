function b = iszero(x)
% ISZERO Check value x.
%
%   b = ISZERO(x) returns whether x is empty or zero.
%
%   See also isempty.

%  Copyright (c) 2016-2017 Yingzhou Li and Chao Yang,
%                          Stanford University and Lawrence Berkeley
%                          National Laboratory
%  This file is distributed under the terms of the MIT License.

b = isempty(x) || (sumel(x ~= 0) == 0);

end