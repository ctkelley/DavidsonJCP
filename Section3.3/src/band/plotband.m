function ax = plotband(cry,rho,efermi,nbnd,varargin)
% PLOTBAND plots the band structure.
%   plotband(cry,rho,efermi,nbnd,endpts) plot the band structure of cry and
%   rho. efermi denotes the fermi level and is optional. If efermi is
%   left empty, the fermi level would not be plotted. nbnd is the number of
%   bands in the plot. endpts is a nkpt by 4 cell. Each row of endpts
%   denotes a three dimensional end points together with a point name, e.g,
%       endpts = { 0   0   0   '\Gamma'
%                  0   0.5 0.5 'X' };
%
%   plotband(cry,rho,efermi,nbnd,endpts,namepts) plot the band structure of
%   cry and rho. endpts is a nkpt by 3 matrix with each row being a three
%   dimensional end points and the corresponding name given in namepts.
%
%   plotband(...,ylimit,gap,hlineflag) ylimit specify the plot range in y
%   axis. For gap < 1, gap denotes the gap between contigous k-points in a
%   line; for gap >= 1, gap denotes the number of point in a line; for
%   vector gap, the first nendkpts-1 values are used as scalar gap before.
%   hlineflag indicates whether plot the vertical lines for each end point.
%
%   See also eband.

%  Copyright (c) 2016-2017 Yingzhou Li, Lin Lin and Chao Yang,
%                          Stanford University,
%                          University of California, Berkeley
%                          and Lawrence Berkeley National Laboratory
%  This file is distributed under the terms of the MIT License.



if iscell(varargin{1})
    nendpts = size(varargin{1},1);
    endpts = cell2mat(varargin{1}(:,1:3));
    namepts = varargin{1}(:,4);
    offset = 1;
else
    nendpts = size(varargin{1},1);
    endpts = varargin{1};
    namepts = varargin{2};
    offset = 2;
end

ylimit = [];
gap = 0.01;
hlineflag = true;

if nargin - 4 - offset > 0
    ylimit = varargin{offset+1};
end

if nargin - 4 - offset > 1
    gap = varargin{offset+2};
end

if nargin - 4 - offset > 2
    hlineflag = varargin{offset+3};
end

if max(size(gap)) == 1
    gap = gap*ones(nendpts,1);
end

npt = zeros(1,nendpts);
ebands = cell(nendpts-1,1);
for it = 1:nendpts-1
    stpt = endpts(it,:);
    edpt = endpts(it+1,:);
    len = sqrt(sum((edpt-stpt).^2));
    if gap(it) > 1
        kpts = linspace(0,len,gap(it))'*(edpt-stpt)/len;
    else
        kpts = (0:gap(it):len)'*(edpt-stpt)/len;
    end
    kpts = kpts + repmat(stpt,size(kpts,1),1);
    npt(it+1) = size(kpts,1) + npt(it);
    ebands{it} = eband(cry,rho,nbnd,kpts);
end
ebands = cell2mat(ebands);

ax = axes;
hold on;
xline = linspace(0,1,size(ebands,1));
for itp = 1:nbnd
    plot(xline,ebands(:,itp),'b-');
end

if isempty(ylimit)
    ylimit = [minel(ebands) maxel(ebands)];
end


xtickval = zeros(nendpts,1);
for it = 1:nendpts
    xtickval(it) = npt(it)/npt(end);
end
set(ax,'XTick',xtickval);
set(ax,'XTickLabel',namepts);

if ~isempty(efermi)
    plot([0,1],[efermi efermi],'k.-');
end

if hlineflag
    for it = 1:nendpts
        plot([npt(it)/npt(end),npt(it)/npt(end)],ylimit,'k-');
    end
end

ylim(ylimit);
ylabel('Energy/Ha')

end