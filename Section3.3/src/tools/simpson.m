function intf = simpson(nMesh,func,rab,dim)
% SIMPSON  Compute the integration of the given function using Simpson's
%          rule.
%    intf = SIMPSON(nMesh,func,rab) integrates function func on mesh. nMesh
%    is the number of mesh points. func provides the function values on the
%    mesh, func(i) = f(r(i)). rab is the derivative of the mesh, i.e.,
%    rab(i) = r(i) * ( D r(i) / D i ) * D i. The output intf is the
%    integration of func, intf = \sum_i c_i func(i)*rab(i) = \int_0^\infty
%    f(r) D r, where c_i is simpson coefficient. The explicit formula is
%    intf = 1/3 * \sum_{j=1}^{j=nMesh/2} [func(2j-2)rab(2j-2) +
%           4 * func(2j-1)rab(2j-1) + func(2j)rab(2j)].
%
%    intf = SIMPSON(nMesh,func,rab,dim) integrates function func on mesh
%    along dim direction. Other input and output parameters are the same
%    as previous one.
%
%    Remark: Currently, nMesh must be a odd number which agrees with
%    QE-5.4.0. If any even number is given, the last point is ignored.
% 
%    See also xmlread, upf.
%
%    Reference page in Help browser
%       doc upfread

%  Copyright (c) 2015-2016 Yingzhou Li and Chao Yang,
%                          Stanford University and Lawrence Berkeley
%                          National Laboratory
%  This file is distributed under the terms of the MIT License.

if nargin < 4
    dim = 0;
elseif dim > ndims(func)
    intf = func;
    return;
else
    dim = dim - 1;
end

ndim = ndims(func);

func = shiftdim(func,dim);

siz = size(func);
ldim = siz(1);
rab = reshape(rab,ldim,1);

idx = 1:2:nMesh-2;
f = func.*repmat(rab,1,siz(2:end));
intf = sum(1/3 * ( f(idx,:) + 4*f(idx+1,:) + f(idx+2,:) ),1);

intf = shiftdim(intf,ndim-dim);

end
