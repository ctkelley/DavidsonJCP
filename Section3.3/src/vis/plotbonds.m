function hfig = plotbonds(varargin)

RC = 0.2; % radius of bonds
FACTOR_CONNECTIVITY = 1.05;

if isa(varargin{1},'Molecule')
    hfig = axes;
    hfig.Color = [0 0 0]; % set background to be black
    mol = varargin{1};
else
    hfig = varargin{1};
    mol = varargin{2};
end

anums = arrayfun(@(it)mol.atoms(it).anum,mol.alist);
xyzlist = mol.xyzlist;
Natoms = size(xyzlist,1);

rads = zeros(Natoms,1);
for it = 1:Natoms
    rads(it) = atomradius(anums(it));
end

% all combinations of atom pairs
pmat = repmat(1:Natoms,Natoms,1);
tpmat = pmat';
pairs = [pmat(:) tpmat(:)];
idx = pairs(:,1)>pairs(:,2);
pairs = pairs(idx,:);
% all interatomic distances
ds = sqrt(sum((xyzlist(pairs(:,1),:) - xyzlist(pairs(:,2),:)).^2,2));
cutds = (rads(pairs(:,1))+rads(pairs(:,2)))*FACTOR_CONNECTIVITY;

% find bonds based on distances of atoms
ks = find(ds<cutds)';

% draw cylinders for each bond
for k = ks % draw sticks for all bounds
    r1 = xyzlist(pairs(k,1),:); % coordinates atom 1
    r2 = xyzlist(pairs(k,2),:); % coordinates atom 2
    
    % bond angles in spherical coordinates
    v = (r2-r1)/norm(r2-r1);
    phi = atan2(v(2),v(1));
    theta = -asin(v(3));
    rotphi = [ cos(phi), -sin(phi), 0;
               sin(phi),  cos(phi), 0;
                      0,         0, 1 ];
    rottheta = [  cos(theta), 0, sin(theta);
                           0, 1,          0;
                 -sin(theta), 0, cos(theta) ];
    rotmat = rotphi*rottheta;
                  
    
    % bond distance minus sphere radii
    bd = ds(k);
    cyl1 = bd/2; % length full bond cylinder
    cyl2 = bd/2; % length half bond cylinder
    
    
    % prototype cylinders for bond
    [z1,y1,x1] = cylinder(RC); % full bond cylinder
    x1(2,:) = x1(2,:) * cyl1; % adjust length
    [z2,y2,x2] = cylinder(RC); % half bond cylinder, thicker
    x2(2,:) = x2(2,:) * cyl2; % adjust length
    
    % rotate cylinders to match bond vector v
    for kk = 1:numel(x1)
        vr = [x1(kk); y1(kk); z1(kk);];
        vr = rotmat*vr;
        x1(kk) = vr(1);
        y1(kk) = vr(2);
        z1(kk) = vr(3);
        
        vr = [x2(kk); y2(kk); z2(kk);];
        vr = rotmat*vr;
        x2(kk) = vr(1);
        y2(kk) = vr(2);
        z2(kk) = vr(3);
    end
    
    % get colors of both atoms
    col1 = atomcolor(anums(pairs(k,1)));
    col2 = atomcolor(anums(pairs(k,2)));
    
    % bond color 1
    surface('XData',r1(1)+x1, ...
        'YData',r1(2)+y1, ...
        'ZData',r1(3)+z1, ...
        'Parent', hfig, ...
        'FaceColor',col1, 'EdgeColor','none', ...
        'UserData', 'bond' );
    
    % bond color 2
    surface('XData',r2(1)-x2, ...
        'YData',r2(2)-y2, ...
        'ZData',r2(3)-z2, ...
        'Parent', hfig, ...
        'FaceColor',col2, 'EdgeColor','none', ...
        'UserData', 'bond' );
end

end