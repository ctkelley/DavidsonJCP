function plotall(varargin)

if isa(varargin{1},'Molecule')
    hfig = axes;
    hfig.Color = [0 0 0]; % set background to be black
    off = 0;
else
    hfig = varargin{1};
    off = 1;
end

mol = varargin{1+off};
H = varargin{2+off};
isoval = varargin{3+off};
alphaval = varargin{4+off};
fdensity = varargin{5+off};
fcenterdensity = varargin{6+off};

plotsetup(hfig,mol);
plotbonds(hfig,mol);
plotatoms(hfig,mol,'bond');
plotatoms(hfig,mol);

if fdensity
    plotdensity(hfig,mol,H,isoval,alphaval,fcenterdensity);
end
axis equal

end
