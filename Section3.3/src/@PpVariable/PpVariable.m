classdef PpVariable
    % PPVARIABLE  Class for transformed pseudopotential data in G-space.
    %    var = PPVARIABLE(mol) reads the pseudopotential provided by users
    %    into PPDATA and transforms the data in real space to G-space.
    %
    %    The PPVARIABLE class contains the following fields.
    %        Field       Explaination
    %      ----------------------------------------------------------
    %        vlocg    	 Local pseudopotential in G-space
    %        rhog        Density read from PpData in G-space
    %        vnltab      Non-local pseudopotential interpolation table
    %        vnlD        Square D matrix for non-local pseudopotential
    %        lll         The degree in the non-local pseudopotential
    %        venums      The numbers of valence electrons
    %        funct       The functionals of the PpData
    %
    %    See also PpData.
    
    %  Copyright (c) 2016-2017 Yingzhou Li and Chao Yang,
    %                          Stanford University and Lawrence Berkeley
    %                          National Laboratory
    %  This file is distributed under the terms of the MIT License.
    properties
        vlocg
        rhog
        vnltab
        vnlD
        lll
        venums
        funct
    end
    methods
        function var = PpVariable(mol)
            
            atoms = mol.atoms;
            ntypes = length(atoms);
            
            % Specify ecut and ecut2, which should be combined in the
            % future
            ecut2    = mol.ecut2;
            if isa(mol,'Crystal')
                ecut = (sqrt(mol.ecut) ...
                    +sqrt(max(sum(mol.kpts.^2,2))/(2*meDef())))^2;
            else
                ecut = mol.ecut;
            end
            
            grid2 = Ggrid(mol,ecut2);
            ng2    = grid2.ng;
            gkk2   = grid2.gkk;
            
            var.vlocg = zeros(ng2,ntypes);
            var.rhog = zeros(ng2,ntypes);
            var.vnltab = cell(1,ntypes);
            var.vnlD = cell(1,ntypes);
            var.lll = cell(1,ntypes);
            var.venums = zeros(1,ntypes);
            
            for it = 1:ntypes
                pp = PpData(atoms(it));
                [var.vlocg(:,it),var.rhog(:,it)] = vloc2g(pp,gkk2);
                [var.vnltab{it},var.vnlD{it}] = vnl2g(pp,ecut);
                var.lll{it} = pp.nonloc.lll;
                var.venums(it) = pp.venum;
                if isempty(var.funct)
                    var.funct = pp.info.functional;
                elseif ~strcmpi(var.funct,pp.info.functional)
                    warning(['The functionals in pseudo potential' ...
                        ' files are not consist, ' var.funct ...
                        ' is used anyway']);
                end
            end
            
        end
    end
end
