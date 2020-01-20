classdef Atom
    % ATOM KSSOLV class for atom
    %    a = ATOM(asym) returns a atom class of the given atom symbol.
    %
    %    a = ATOM(anum) returns a atom class of the given atom number.
    %
    %    The atom class contains the following fields.
    %        Field     Explaination
    %      ----------------------------------------------------------
    %        symbol    Atom symbol
    %        anum      Atom number
    %        amass     Atom mass
    %        venum     Number of valence electrons
    %        iloc      TBA
    %        occs      TBA
    %        occp      TBA
    %        occd      TBA
    %        iso       TBA
    %        ic        TBA
    %        isref     TBA
    %        ipref     TBA
    %        idref     TBA
    %
    %    See also Molecule.
    
    %  Copyright (c) 2016-2017 Yingzhou Li and Chao Yang,
    %                          Stanford University and Lawrence Berkeley
    %                          National Laboratory
    %  This file is distributed under the terms of the MIT License.
    properties (SetAccess = protected)
        symbol
        anum
        amass
        venum
        iloc
        occs
        occp
        occd
        iso
        ic
        isref
        ipref
        idref
    end
    methods
        function a = Atom(asym)
            if nargin < 1
                return;
            end
            if ( ischar(asym) )
                num = sym2num(asym);
            else
                num = asym;
            end
            a.symbol = num2sym(num);
            a.anum = num;
            [a.amass,a.venum,a.iloc,a.occs,a.occp,a.occd,a.iso,a.ic,...
                a.isref,a.ipref,a.idref] = atominfo(num);
        end
    end
    methods (Static)
        function z = zeros(varargin)
            if (nargin == 0)
                z = Atom;
            elseif any([varargin{:}] <= 0)
                z = Atom.empty(varargin{:});
            else
                z = repmat(Atom,varargin{:});
            end
        end
    end
end