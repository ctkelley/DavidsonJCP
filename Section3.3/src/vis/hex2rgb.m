function rgb = hex2rgb(hex,range)
% HEX2RGB Converts hex color values to RGB arrays.
%   rgb = HEX2RGB(hex) returns RGB color values in an n x 3 array given the
%   n strings in hex. Values are scaled from 0 to 1 by default.
%
%   rgb = HEX2RGB(hex,range) returns RGB color values in an n x 3 array
%   given the n strings in hex. Values are scaled from 0 to range-1.
%   Usually, range=255 is specified.
%
%   Examples:
%    >> rgb = hex2rgb('#334D66')
%    rgb =
%        0.2000    0.3020    0.4000
%
%
% 	 >> rgb = hex2rgb('334D66')  % <-the # sign is optional
%    rgb =
%        0.2000    0.3020    0.4000
%
%
%    >> hexs = ['#334D66';'#8099B3'];
% 	 >> rgbs = hex2rgb(hexs)
%    rgbs =
%        0.2000    0.3020    0.4000
%        0.5020    0.6000    0.7020
%
%    >> hexs = {'#334D66';'#8099B3'};
% 	 >> rgbs = hex2rgb(hexs)
%    rgbs =
%        0.2000    0.3020    0.4000
%        0.5020    0.6000    0.7020
%
%   See also rgb2hex, dec2hex, hex2num, and ColorSpec.

%  Copyright (c) 2016-2017 Yingzhou Li and Chao Yang,
%                          Stanford University and Lawrence Berkeley
%                          National Laboratory
%  Copyright (c) 2014 Chad A. Greene
%  This file is distributed under the terms of the MIT License.

if iscell(hex)
    if isrow(hex)
        hex = hex';
    end
    % If input is cell, convert to matrix:
    hex = cell2mat(hex);
end

if strcmpi(hex(1,1),'#')
    hex(:,1) = [];
end

if nargin == 1
    range = 1;
end

% Convert from hex to rgb:
if range == 1
    rgb = reshape(sscanf(hex.','%2x'),3,[]).'/255;
else
    rgb = round(reshape(sscanf(hex.','%2x'),3,[]).'*(range/255));
end

end