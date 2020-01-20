function hfig = plotcell(varargin)

if isa(varargin{1},'Molecule')
    hfig = axes;
    hfig.Color = [0 0 0]; % set background to be black
    mol = varargin{1};
else
    hfig = varargin{1};
    mol = varargin{2};
end

C = mol.supercell;
vert = nan(8,3);
fac = nan(6,4);

vert(1,:) = [ 0 0 0];
vert(2,:) = C(1,:);
vert(3,:) = C(2,:);
vert(4,:) = C(3,:);
vert(5,:) = C(1,:) + C(2,:);
vert(6,:) = C(1,:) + C(3,:);
vert(7,:) = C(2,:) + C(3,:);
vert(8,:) = C(1,:) + C(2,:) + C(3,:);

fac(1,:) = [1 2 5 3];
fac(2,:) = [1 2 6 4];
fac(3,:) = [1 3 7 4];
fac(4,:) = [2 5 8 6];
fac(5,:) = [3 5 8 7];
fac(6,:) = [4 6 8 7];

patch('Vertices',vert, 'Faces',fac, ...
    'EdgeColor',[1 1 1], 'LineWidth',2, ...
    'UserData','cell');

end