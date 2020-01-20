function H = subsref(BH,S)
% BLOCHHAM/SUBSREF Subsref function for BlochHam class
%    H = BH{ik} returns the Ham class with ik-th k-point.
%
%    See also Ham, BlochHam.

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
        H = Ham();
        H.n1      = BH.n1;
        H.n2      = BH.n2;
        H.n3      = BH.n3;
        H.gkin    = BH.gkincell{ik};
        H.idxnz   = BH.idxnz;
        H.vtot    = BH.vtot;
        H.vion    = BH.vion;
        H.vext    = BH.vext;
        H.vnp     = BH.vnp;
        H.vnlmat  = BH.vnlmatcell{ik};
        H.vnlsign = BH.vnlsigncell{ik};
        H.rho     = BH.rho;
    otherwise
        H = builtin('subsref',BH,S);
end

end