function [idxev,ev] = sortev(ev,n)
% SORTEV returns the smallest n eigenvalues in ev.
%   [idxev,ev] = SORTEV(ev,n) find the smallest n eigenvalues in ev and
%   returns the corresponding indices and eigenvalues.
%
%   See also sort.

%  Copyright (c) 2016-2017 Yingzhou Li and Chao Yang,
%                          Stanford University and Lawrence Berkeley
%                          National Laboratory
%  This file is distributed under the terms of the MIT License.

nkpts = numel(ev);
allev = cell2mat(ev);
allev = sort(allev);
if n < numel(allev)
    evcut = (allev(n)+allev(n+1))/2;
else
    evcut = allev(end)+1;
end

idxev = cell(nkpts,1);
for ik = 1:nkpts
    idxev{ik} = find(ev{ik} < evcut);
    ev{ik} = ev{ik}(idxev{ik});
end

end