function jl = bessel(l,x)
% BESSEL  Compute the approximated Bessel function of the first kind.
%    jl = BESSEL(l,x) calculates the approximated Bessel function of the
%    first kind. It is suppose to implement the similar function as
%    besselj in Matlab. However, the output is very different from the
%    Matlab default bessel functions. The goal of this implementation is to
%    match the sph_bes implementation in Quantum Espresso. We currently
%    only support the positive order l.
% 
%    See also besseli, besselj, besselk, besselh.

%  Copyright (c) 2015-2016 Yingzhou Li and Chao Yang,
%                          Stanford University and Lawrence Berkeley
%                          National Laboratory
%  This file is distributed under the terms of the MIT License.

xseries = 0.05;

jl = zeros(size(x));

idxzero = abs(x)<eps;
if l == 0
    jl(idxzero) = 1;
else
    jl(idxzero) = 0;
end

idx = (abs(x)>=eps) & (abs(x) <= xseries);
xcur = x(idx);
xlcur = xcur.^l;
jl(idx) = xlcur/doublefactorial(2*l+1) ...
        .* ( 1 - xcur.^2/1/2/(2*l+3) ) ...
        .* ( 1 - xcur.^2/2/2/(2*l+5) ) ...
        .* ( 1 - xcur.^2/3/2/(2*l+7) ) ...
        .* ( 1 - xcur.^2/4/2/(2*l+9) );

idxseries = abs(x) > xseries;
xcur = x(idxseries);

switch l
    case 0
        jl(idxseries) = sin(xcur) ./ xcur;
    case 1
        jl(idxseries) = ( sin(xcur) ./ xcur - cos(xcur) )./ xcur;
    case 2
        jl(idxseries) = ( sin(xcur) .* ( 3./xcur - xcur ) ...
                        - cos(xcur) .* 3 ...
                        ) ./ xcur.^2;
    case 3
        jl(idxseries) = ( sin(xcur) .* ( 15./xcur - 6 * xcur ) ...
                        + cos(xcur) .* ( xcur.^2 - 15 ) ...
                        ) ./ xcur.^3;
    case 4
        jl(idxseries) = ( sin(xcur) .* ( 105 - 45*xcur.^2 + xcur.^4 ) ...
                        + cos(xcur) .* ( 10*xcur.^3 - 105*cur ) ...
                        ) ./ xcur.^5;
    case 5
        sx = sin(xcur);
        cx = cos(xcur);
        jl(idxseries) = ( -cx ...
                      - ( 945*cx./xcur.^4 + 105*cx./xcur.^2 ...
                        + 945*sx./xcur.^5 - 420*sx./xcur.^3 ...
                        +  15*sx./xcur ) ...
                        ) ./ xcur;
    case 6
        sx = sin(xcur);
        cx = cos(xcur);
        jl(idxseries) = ( -10395*cx./xcur.^5 +  1260*cx./xcur.^3 ...
                          -   21*cx./xcur ...
                          -      sx          + 10395*sx./xcur.^6 ...
                          - 4725*sx./xcur.^4 +   210*sx./xcur.^2 ...
                        ) ./ xcur;
    otherwise
        error('Bessel function does not support the given l.')
end

end