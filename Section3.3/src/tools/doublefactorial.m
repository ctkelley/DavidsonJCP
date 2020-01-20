function k = doublefactorial(n)
% DOUBLEFACTORIAL double facotorial function.
%   k = DOUBLEFACTORIAL(n) for integer n, is the double factorial of n,
%   i.e., k = n*(n-2)*(n-4)*...*2 for even n, and k = n*(n-2)*...*1 for odd
%   n.
%
%   See also factorial.

%  Copyright (c) 2015-2016 Yingzhou Li and Chao Yang,
%                          Stanford University and Lawrence Berkeley
%                          National Laboratory
%  This file is distributed under the terms of the MIT License.

k = prod(n:-2:1);

end