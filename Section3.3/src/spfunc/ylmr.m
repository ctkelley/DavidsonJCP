function YlmrTab = ylmr(x,y,z,maxl)
% YLMR creates a real spherical harmonics table.
%
%   YlmrTab = YLMR(x,y,z,maxl) returns the real spherical harmonics table
%   for l in the range of [0,1,...,maxl] at 3 dimensional locations given
%   by the column vectors x, y, z. The return table is of size n by
%   (1+maxl)^2. Here n is the number of points. Each row of the return
%   table corresponds to a location point. Meanwhile, the columns indicate
%   the degree l and order m, (l,m), in a natural increaseing order, i.e.,
%   (0,0), (1,-1), (1,0), (1,1), (2,-2), ..., (maxl,maxl-1), (maxl,maxl).
%   It is suggested that maxl should be smaller than 8, otherwise, the
%   numerical results are not stable.
%
%   YlmrTab = YLMR(x,y,z) returns the real spherical harmonics table for l
%   in the range of [0,1,...,9] at 3 dimensional locations given by the
%   column vectors x, y, z. The return table is of size n by 100 of the
%   same structure as previous description.
%
%   YlmrTab = YLMR(xyz,maxl) returns the real spherical harmonics table for
%   l in the range of [0,1,...,maxl] at 3 dimensional locations given by
%   the column vectors xyz, which is of dimension n by 3. Here n is the
%   number of points. The return table is the same as previous description.
%
%   YlmrTab = YLMR(xyz) returns the real spherical harmonics table for l
%   in the range of [0,1,...,9] at 3 dimensional locations given by the
%   column vectors xyz, which is of dimension n by 3. Here n is the number
%   of points. The return table is the same as previous description.
%
%   See also ylm.

%  Copyright (c) 2016-2017 Yingzhou Li and Chao Yang,
%                          Stanford University and Lawrence Berkeley
%                          National Laboratory
%  This file is distributed under the terms of the MIT License.

switch nargin
    case 1,
        y = x(:,2);
        z = x(:,3);
        x = x(:,1);
        maxl = 9;
    case 2,
        y = x(:,2);
        z = x(:,3);
        x = x(:,1);
    case 3,
        maxl = 9;
end

len = length(x);

if maxl == 0
    YlmrTab = repmat(sqrt(1/4/pi),len,1);
    return;
end

r = sqrt(x.^2+y.^2+z.^2);
zeroidxr = abs(r)<eps;
zeroidxx = abs(x)<eps;
cost = z./r;
cost(zeroidxr) = 0;

phi = atan(y./x);
lzidx = x<-eps;
phi(zeroidxx) = pi/2*sign(y(zeroidxx));
phi(lzidx) = phi(lzidx)+pi;

sent = sqrt(max(0,1-cost.^2));

YlmrTab = zeros(len,(1+maxl)^2);
Q = zeros(len,maxl+1,maxl+1);

lm = 0;
for ll = 0:maxl
    
    pref = sqrt((2*ll+1)/4/pi);
    
    if ll == 0
        Q(:,ll+1,1) = ones(len,1);
    elseif ll == 1
        Q(:,ll+1,1) = cost;
        Q(:,ll+1,2) = -sent/sqrt(2);
    else
        for m = 0:ll-2
            Q(:,ll+1,m+1) = cost*(2*ll-1)/sqrt(ll^2-m^2).*Q(:,ll,m+1) ...
                - sqrt((ll-1)^2-m^2)/sqrt(ll^2-m^2)*Q(:,ll-1,m+1);
        end
        
        Q(:,ll+1,ll) = cost*sqrt(2*ll-1).*Q(:,ll,ll);
        Q(:,ll+1,ll+1) = -sqrt(2*ll-1)/sqrt(2*ll)*sent.*Q(:,ll,ll);
    end
    
    % m < 0
    for m = -ll:-1
        lm = lm+1;
        YlmrTab(:,lm) = pref*sqrt(2)*Q(:,ll+1,-m+1).*sin(-m*phi);
    end
    
    % m = 0
    lm = lm+1;
    YlmrTab(:,lm) = pref*Q(:,ll+1,1);
   
    % m > 0
    for m = 1:ll
        lm = lm+1;
        YlmrTab(:,lm) = pref*sqrt(2)*Q(:,ll+1,m+1).*cos(m*phi);
    end
end

end
