function hfig = plotdensity(varargin)
%ISORHO Plot isosurface of charge density
%   h = ISORHO(rho,isoval) plots the isosurface of the electron density rho
%   at isoval and returns the graphics handle h. rho is first FFTshifted to
%   the center.
%
%   h = ISORHO(rho,isoval,'noshift') overrides the shift behavior
%
%   See also: ISOPSI

if isa(varargin{1},'Molecule')
    hfig = axes;
    off = 0;
else
    hfig = varargin{1};
    off = 1;
end

mol = varargin{off+1};
H = varargin{off+2};
isoval = varargin{off+3};

if nargin-off == 3
    alphaval = 1;
    centered = true;
else
    alphaval = varargin{off+4};
    if nargin-off == 4
        centered = true;
    else
        centered = varargin{off+5};
    end
end

n1 = mol.n1;
n2 = mol.n2;
n3 = mol.n3;

if centered
    if ~isa(H,'Ham') && ~isa(H,'BlochHam')
        rho = fftshift(H);
    else
        rho = fftshift(H.rho);
    end
    [X,Y,Z] = ndgrid(floor(-n1/2):(n1-2)/2, ...
        floor(-n2/2):(n2-2)/2, ...
        floor(-n3/2):(n3-2)/2 );
else
    if ~isa(H,'Ham') && ~isa(H,'BlochHam')
        rho = H;
    else
        rho = H.rho;
    end
    [X,Y,Z] = ndgrid(0:n1-1,0:n2-1,0:n3-1);
end

XYZ = [X(:)/n1 Y(:)/n2 Z(:)/n3]*mol.supercell;
X = reshape(XYZ(:,1),n1,n2,n3);
Y = reshape(XYZ(:,2),n1,n2,n3);
Z = reshape(XYZ(:,3),n1,n2,n3);

patch(isosurface(X,Y,Z,rho,isoval), ...
    'edgecolor', 'none', 'Parent',hfig, 'facecolor', 'g', ...
    'UserData', 'density', 'FaceAlpha', alphaval );

end