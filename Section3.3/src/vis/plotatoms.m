function hfig = plotatoms(varargin)

if isa(varargin{1},'Molecule')
    hfig = axes;
    hfig.Color = [0 0 0]; % set background to be black
    mol = varargin{1};
    if nargin == 2
        radtype = varargin{2};
    else
        radtype = 'unified';
    end
else
    hfig = varargin{1};
    mol = varargin{2};
    if nargin == 3
        radtype = varargin{3};
    else
        radtype = 'unified';
    end
end

if isa(mol,'Atom')
    % TODO
else
    anums = arrayfun(@(it)mol.atoms(it).anum,mol.alist);
    if strcmpi(radtype,'unified')
        aradius = 0.5*ones(size(mol.alist));
    elseif strcmpi(radtype,'bond')
        aradius = 0.2*ones(size(mol.alist));
    else
        aradius = arrayfun(@(it)atomradius(mol.atoms(it).anum),mol.alist);
    end
    xyzlist = mol.xyzlist;
    C = mol.supercell;
end

% basic sphere
[sx,sy,sz] = sphere;

for it = 1:size(xyzlist,1)
    col = atomcolor(anums(it));
    rad = aradius(it);
    % draw sphere
    if strcmpi(radtype,'bond')
        surface('XData', xyzlist(it,1)+rad*sx, ...
            'YData', xyzlist(it,2)+rad*sy, ...
            'ZData', xyzlist(it,3)+rad*sz, ...
            'Parent', hfig, ...
            'FaceColor', col , 'EdgeColor', 'none', ...
            'UserData', 'bond' );
    else
        surface('XData', xyzlist(it,1)+rad*sx, ...
            'YData', xyzlist(it,2)+rad*sy, ...
            'ZData', xyzlist(it,3)+rad*sz, ...
            'Parent', hfig, ...
            'FaceColor', col , 'EdgeColor', 'none', ...
            'UserData', 'atom' );
    end
end

end