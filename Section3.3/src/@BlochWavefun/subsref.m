function X = subsref(BX,S)
% BLOCHWAVEFUN/SUBSREF Subsref function for Bloch wave function class
%    X = BX{ik} returns the a single wave function.
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
            X = builtin('subsref',BX.wavefuncell{ik},S(2:end));
        else
            X = BX.wavefuncell{ik};
        end
    otherwise
        X = builtin('subsref',BX,S);
end

end