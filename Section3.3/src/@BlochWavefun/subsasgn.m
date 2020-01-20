function BX = subsasgn(BX,S,X)
% BLOCHWAVEFUN/SUBSASGN Subsasgn function for Bloch wave function class
%    BX{ik} = X assign the single wave function.
%
%    See also Wavefun, BlochWavefun.

%  Copyright (c) 2016-2017 Yingzhou Li and Chao Yang,
%                          Stanford University and Lawrence Berkeley
%                          National Laboratory
%  This file is distributed under the terms of the MIT License.

switch S(1).type
    case '{}'
        ik = S(1).subs{1};
        if numel(ik) > 1
            error('Wrong sub index.')
        end
        if numel(S) > 1
            BX.wavefuncell{ik} = ...
                builtin('subsasgn',BX.wavefuncell{ik},S(2:end),X);
        else
            BX.wavefuncell{ik} = X;
        end
    otherwise
        BX = builtin('subsasgn',BX,S,X);
end

end