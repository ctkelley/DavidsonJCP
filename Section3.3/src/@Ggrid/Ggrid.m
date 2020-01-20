classdef Ggrid
    % GGRID is a class that defines the reciprocal space grid. The properties in the
    %   class are all protected, i.e., direct editing is prohibited. The
    %   properties are explained in the following table.
    %     ============================================================
    %     Property     Comment
    %     -----------  -----------------------------------------------
    %     ecut         Energy cut.
    %     ng           Number of the valid points.
    %     idxnz        List of valid indices w.r.t the original index
    %                  list.
    %     gkk          An array of the distance square of the valid
    %                  points in G-space.
    %     gkx,gky,gkz  Arraies of (x,y,z)-value of the valid points in
    %                  G-space
    %     ============================================================
    %
    %   Grid = GGRID returns an empty frequency mask.
    %
    %   Grid = GGRID(mol) returns a frequency mask for the mol, which
    %   is an instance of Molecule class in KSSOLV. The mask is associated
    %   with the mol with the same grid size and energy cut. Shift is not
    %   added.
    %
    %   Grid = GGRID(mol,ecut) returns a frequency mask for the mol
    %   with the same grid size but the given energy cut, ecut. Shift is
    %   not added.
    %
    %   Grid = GGRID(gx,gy,gz,ecut) returns a frequency mask for
    %   locations (gx,gy,gz) with the given energy cut. Shift is not added.

    %  Copyright (c) 2016-2017 Yingzhou Li and Chao Yang,
    %                          Stanford University and Lawrence Berkeley
    %                          National Laboratory
    %  This file is distributed under the terms of the MIT License.
    
    properties (SetAccess = protected)
        ecut
        ng
        idxnz
        gkk
        gkx
        gky
        gkz
    end
    methods
        function Grid = Ggrid(varargin)
            
            %===========================================================
            % Input process
            
            % Return empty mask
            if nargin == 0
                return;
            end
            
            if isa(varargin{1},'Molecule')
                % Process input with mol
                mol = varargin{1};
                n1 = mol.n1;
                n2 = mol.n2;
                n3 = mol.n3;
                [I,J,K] = ndgrid((0:n1-1)-((0:n1-1) >= n1/2)*n1, ...
                    (0:n2-1)-((0:n2-1) >= n2/2)*n2, ...
                    (0:n3-1)-((0:n3-1) >= n3/2)*n3);
                pregkx = I(:);
                pregky = J(:);
                pregkz = K(:);
                C  = mol.supercell;
                
                if nargin > 1 && isnumeric(varargin{2})
                    Gecut = varargin{2};
                else
                    Gecut = mol.ecut;
                end
                
            elseif isnumeric(varargin{1}) && isnumeric(varargin{2}) ...
                    && isnumeric(varargin{3}) && isnumeric(varargin{4})
                % Process input with locations
                pregkx = varargin{1};
                pregky = varargin{2};
                pregkz = varargin{3};
                C  = eye(3);
                
                Gecut = varargin{4};
                
            else
                error('Invalid input.');
            end
            % End of input process
            
            %===========================================================
            % Construction of Ggrid from Molecule
            %
            Grid.ecut = Gecut;
            xyz = [pregkx pregky pregkz];
            
            % Define the nonshifted mask
            % NOTE: the transpose is very important.
            Creci = 2*pi*inv(C)';

            kkxyz = xyz*Creci;
            kk = sum(kkxyz.*kkxyz,2);

            Grid.idxnz = find(kk<=Grid.ecut*2*meDef());
            
            Grid.ng    = length(Grid.idxnz);
            Grid.gkk   = kk(Grid.idxnz);
            Grid.gkx   = kkxyz(Grid.idxnz,1);
            Grid.gky   = kkxyz(Grid.idxnz,2);
            Grid.gkz   = kkxyz(Grid.idxnz,3);
        end
    end
end
