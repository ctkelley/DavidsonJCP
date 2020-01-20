function [pathfile,file] = dirregexp(path,exp)
% DIRREGEXP returns the first file under the path that matches the regular
%   expression.
%   [pathfile,file] = DIRREGEXP(path,exp) find all files under the path,
%   and return the first one that matches the regular expression exp. If
%   none of the files matches the exp, the function will return an empty
%   string. pathfile is the file with path whereas file is the single file
%   name.
%
%   See also regexp, dir.

%  Copyright (c) 2016-2017 Yingzhou Li and Chao Yang,
%                          Stanford University and Lawrence Berkeley
%                          National Laboratory
%  This file is distributed under the terms of the MIT License.

dirlist = dir(path);
for it = 1:length(dirlist)
    matchStr = regexp(dirlist(it).name,exp,'match');
    if ~isempty(matchStr)
        file = matchStr{1};
        pathfile = [path matchStr{1}];
        return;
    end
end
file = '';
pathfile = '';

end