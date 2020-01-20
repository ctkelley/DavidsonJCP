function H = updateV(H,mol,rho)
% UPDATEV update of local potential
%   H = UPDATEV(H) updates the total local potential and keeps it band
%   limited. Typically you would do this after changing H.vnp and/or H.vext
%   and tries to compute H*X using an updated H.vtot, for example:
%
%   >> H.vext = ...; % new external potential
%   >> H = updateV(H); % now H.vtot reflects the change in H.vext
%   >> Y = H*X;
%
%   H = UPDATEV(H,mol,rho) first updates H.vnp using the new density rho
%   and then updates H.vtot. Previous changes to H.vnp and H.vext are also
%   incorporated in the update.
%
%   See also: HAM

if nargin == 3
    H.rho = rho;
    % Calculate Hartree & xc potential
    [vhart,vxc,~,~]=getvhxc(mol,abs(rho));
    H.vnp = vhart+vxc;
end
vtot = H.vion + H.vnp + H.vext;
vsum = sumel(vtot) / (H.n1 * H.n2 * H.n3);
% keep vtot band limited (sinc interpolation)
gm2 = zeros(H.n1,H.n2,H.n3);
gm2(H.idxnz2) = 1;
iz = gm2==0;
vfft = fft3(vtot);
vfft(iz) = 0;
vtot = ifft3(vfft);
H.vtot = vtot - vsum;

end