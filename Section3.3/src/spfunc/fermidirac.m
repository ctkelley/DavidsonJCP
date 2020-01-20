function f = fermidirac(ev,efermi,Tbeta)
% FERMIDIRAC calculates Fermi-Dirac distribution.
%    f = FERMIDIRAC(ev,efermi,Tbeta) returns the Fermi-Dirac distribution
%    with Fermi energy efermi, temperature Tbeta at energy ev. The
%    Fermi-Dirac distribution is defined as follows,
%       f(e) = \frac{1}{ exp{ (e-e_{Fermi})/T_{beta} } +1 }.
%   
%    See also Molecule.

%  Copyright (c) 2016-2017 Yingzhou Li and Chao Yang,
%                          Stanford University and Lawrence Berkeley
%                          National Laboratory
%  This file is distributed under the terms of the MIT License.

f = 1./(1+exp((ev-efermi)/Tbeta));

end
