function varargout = atominfo(anum,varargin)
% ATOMINFO  returns the atom information of the given atom.
%   [amass,venum,iloc,occs,occp,occd,iso,ic,isref,ipref,idref] = 
%   ATOMINFO(anum) returns the atom information of the given atom.
%
%   [info1,info2,...] = ATOMINFO(anum,infoname1,infoname2,...) returns the
%   requested atom information. The options for infoname are 'anum',
%   'amass', 'venum', 'iloc', 'occs', 'occp', 'occd', 'iso', 'ic', 'isref',
%   'ipref', 'idref'.
%
%   Example:
%
%       [amass,occs,venum] = atominfo('Li','amass','occs','venum')
%
%       amass =
% 
%           6.9390
% 
% 
%       occs =
% 
%            1
% 
% 
%       venum =
% 
%            1
%
%
%   See also sym2num, num2sym, Atom, GlobalPeriodicTable.

%  Copyright (c) 2016-2017 Yingzhou Li and Chao Yang,
%                          Stanford University and Lawrence Berkeley
%                          National Laboratory
%  This file is distributed under the terms of the MIT License.

if ischar(anum)
   anum = sym2num(anum);
end

if ~isnumeric(anum) || anum <1 || anum > 110
    error('Invalid input');
end

global PeriodicTable;
if isempty(PeriodicTable)
    GlobalPeriodicTable();
end

varargout = cell(1,nargout);
if nargin == 1 && nargout == 11
    for it = 1:11
        varargout{it} = PeriodicTable(anum,it+1);
    end
elseif nargin == nargout+1
    for it = 1:nargout
        switch varargin{it}
            case 'anum'
                varargout{it} = PeriodicTable(anum,1);
            case 'amass'
                varargout{it} = PeriodicTable(anum,2);
            case 'venum'
                varargout{it} = PeriodicTable(anum,3);
            case 'iloc'
                varargout{it} = PeriodicTable(anum,4);
            case 'occs'
                varargout{it} = PeriodicTable(anum,5);
            case 'occp'
                varargout{it} = PeriodicTable(anum,6);
            case 'occd'
                varargout{it} = PeriodicTable(anum,7);
            case 'iso'
                varargout{it} = PeriodicTable(anum,8);
            case 'ic'
                varargout{it} = PeriodicTable(anum,9);
            case 'isref'
                varargout{it} = PeriodicTable(anum,10);
            case 'ipref'
                varargout{it} = PeriodicTable(anum,11);
            case 'idref'
                varargout{it} = PeriodicTable(anum,12);
            otherwise
                error('Invalid input');
        end
    end
else
    error('Invalid input and output');
end

end
