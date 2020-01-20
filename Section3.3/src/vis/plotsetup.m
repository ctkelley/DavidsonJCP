function hfig = plotsetup(varargin)

if nargin == 0
    hfig = axes;
    C = eye(3);
    xyzlist = [0.5 0.5 0.5];
elseif nargin == 1
    if isa(varargin{1},'Molecule')
        hfig = axes;
        C = varargin{1}.supercell;
        xyzlist = varargin{1}.xyzlist;
    else
        hfig = varargin{1};
        C = eye(3);
        xyzlist = [0.5 0.5 0.5];
    end
else
    hfig = varargin{1};
    C = varargin{2}.supercell;
    xyzlist = varargin{2}.xyzlist;
end

xyzlistmin = (fix(xyzlist/C)-1)*C;
xyzlistmax = (fix(xyzlist/C)+1)*C;

xmin = min(xyzlistmin(:,1));
xmax = max(xyzlistmax(:,1));
ymin = min(xyzlistmin(:,2));
ymax = max(xyzlistmax(:,2));
zmin = min(xyzlistmin(:,3));
zmax = max(xyzlistmax(:,3));
vimin = min([xmin,ymin,zmin]);
vimax = max([xmax,ymax,zmax]);

hfig.Color = [0 0 0];

if ~verLessThan('matlab','8.6')
    hfig.XAxis.Visible = 'off';
    hfig.YAxis.Visible = 'off';
    hfig.ZAxis.Visible = 'off';
end

hfig.XLimMode = 'manual';
hfig.YLimMode = 'manual';
hfig.ZLimMode = 'manual';

hfig.XLim = [vimin vimax];
hfig.YLim = [vimin vimax];
hfig.ZLim = [vimin vimax];

hfig.NextPlot = 'add';

hfig.ClippingStyle = 'rectangle';

hfig.View = [-37.5 30];

camlight;

end