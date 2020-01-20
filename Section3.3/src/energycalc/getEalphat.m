function Ealphat = getEalphat(mol)
%
% Calculate the Ealpha energy
%
% usage: Ealphat = getealphat(mol);
%        where mol is a Molecule object;
%

%
% check the validity of the input argument.
%
%...

%fprintf('Calculating Ealphat energy...\n');
%t0      = cputime;
%
% the number of different atom types
%
%
natoms = mol.natoms;
vol    = mol.vol;
ntypes = length(natoms);

ealpha   = zeros(1,ntypes);

for j = 1:ntypes;
    % TODO: the evaluation of the vloc here is not correct.
%    anum = inz(j);
%    pp = PpData(anum);
%    vloc = pp.vloc;
%    r = pp.r;
%    ch = pp.venum;
%    nrr = size(vloc,1);
%    i15 = find(r(2:nrr-1)<15.0)+1;
%    s = sum( (ch*r(i15)+vloc(i15).*r(i15).^2).*(r(i15+1)-r(i15-1))/2 );
%    %
%    ealpha(inz(j)) = s*4*pi;
end
%
Ealphat = sum(ealpha.*natoms);
Ealphat = Ealphat*mol.nel/vol;

end