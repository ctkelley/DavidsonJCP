classdef Wavefun
    % WAVEFUN KSSOLV class for wave function
    %    X = WAVEFUN() returns an empty wave function.
    %
    %    X = WAVEFUN(n1,n2,n3) returns a zero wave function with size n1
    %    by n2 by n3.
    %
    %    X = WAVEFUN(n1,n2,n3,ncols) returns ncols zero wave functions with
    %    size n1 by n2 by n3.
    %
    %    X = WAVEFUN(psi) constructs wave functions by cells of 3D array
    %    psi.
    %
    %    X = WAVEFUN(psi,n1,n2,n3) returns ncols wave functions with
    %    size n1 by n2 by n3 and filled by psi.
    %
    %    X = WAVEFUN(psi,n1,n2,n3,idxnz) constructs compact wave functions
    %    by the psi with size n1 by n2 by n3. And idxnz indicates the
    %    indices for non-zero entries.
    %
    %    The wavefun class contains the following fields.
    %        Field       Explaination
    %      ----------------------------------------------------------
    %        n1,n2,n3    Number of discretization points in each dimension
    %        nrows       Number of rows in the wave function
    %        ncols       Number of columns in the wave function
    %        idxnz       Indices for non-zero entries in the n1,n2,n3
    %        iscompact   Indicator for compact format storage
    %        trans       Indicator for transpose
    %        psi         2D data matrix storing each wave function in each
    %                    column
    %        occ         Occupation rate
    %
    %    See also Atom.
    
    %  Copyright (c) 2016-2017 Yingzhou Li and Chao Yang,
    %                          Stanford University and Lawrence Berkeley
    %                          National Laboratory
    %  This file is distributed under the terms of the MIT License.
    properties (SetAccess = protected)
        n1 = 0
        n2 = 0
        n3 = 0
        idxnz
        iscompact = 0
        trans = 0
        psi
    end
    properties (SetAccess = public)
        occ
    end
    methods
        function X = Wavefun(varargin)
            switch (nargin)
                case 0
                    return;
                case 1
                    Xin = varargin{1};
                    if iscell(Xin)
                        [X.n1,X.n2,X.n3] = size(Xin{1});
                        nrows = numel(Xin{1});
                        ncols = length(Xin);
                        X.psi = zeros(nrows,ncols);
                        for it = 1:length(Xin)
                            X.psi(:,it) = Xin{it}(:);
                        end
                    else
                        if ndims(Xin) == 3
                            [X.n1,X.n2,X.n3]=size(Xin);
                            X.psi  = Xin(:);
                        elseif ndims(Xin) == 4
                            [X.n1,X.n2,X.n3,ncols]=size(Xin);
                            nrows = X.n1*X.n2*X.n3;
                            X.psi  = reshape(Xin,nrows,ncols);
                        end
                    end
                case 3
                    X.n1  = varargin{1};
                    X.n2  = varargin{2};
                    X.n3  = varargin{3};
                case 4
                    if ndims(varargin{1}) > 1
                        X.n1  = varargin{2};
                        X.n2  = varargin{3};
                        X.n3  = varargin{4};
                        X.psi = varargin{1};
                    else
                        X.n1  = varargin{1};
                        X.n2  = varargin{2};
                        X.n3  = varargin{3};
                        nrows  = X.n1*X.n2*X.n3;
                        ncols  = varargin{4};
                        X.psi = zeros(nrows,ncols);
                    end
                case 5
                    Xin = varargin{1};
                    X.n1 = varargin{2};
                    X.n2 = varargin{3};
                    X.n3 = varargin{4};
                    X.idxnz = varargin{5};
                    X.iscompact = 1;
                    if iscell(Xin)
                        X.psi = zeros(numel(Xin{1}),length(Xin));
                        for it = 1:length(Xin)
                            X.psi(:,it) = Xin{it}(:);
                        end
                    else
                        X.psi  = Xin;
                    end
                otherwise
                    error('Wrong number of arguments');
            end
        end
        function ind = end(X,k,n)
            szd = size(X.psi);
            if k < n
                ind = szd(k);
            else
                ind = prod(szd(k:end));
            end
            
        end
    end
end