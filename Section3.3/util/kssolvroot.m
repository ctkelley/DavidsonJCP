function path = kssolvroot()
% KSSOLVROOT Root directory of KSSOLV.
%   S = KSSOLVROOT returns a string that is the absolute path of the
%   directory where the KSSOLVE software is installed.
%  
%   See also matlabroot.
%

%  Copyright (c) 2016-2017 Yingzhou Li and Chao Yang,
%                          Stanford University and Lawrence Berkeley
%                          National Laboratory
%  This file is distributed under the terms of the MIT License.
path = mfilename('fullpath');
path = path(1:end-15);

end