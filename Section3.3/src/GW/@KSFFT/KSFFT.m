function F = KSFFT(varargin)
%
% Fast Fourier transform object for KSSOLV
%
% usage:  F = KSFFT(mol);  constructs an inverse Fourier transform 
%                          that converts Fourier coefficients
%                          of a wavefunction (saved as a vector)
%                          to a wavefunction in real space
%                          saved as a vector
%
%         Once constructed, F can be used as a matrix
%               y = F*x
%         gives the inverse Fourier transform of x 
%
%               x = F'*y;
%         gives the Fourier transform of y
%
%         F contains the volume factor
%

switch (nargin)
  case 0
     F.n1 = [];
     F.n2 = [];
     F.n3 = [];
     F.idxnz = [];
     F.forward = 1;
     F.inverse = 0;
     F = class(F,'KSFFT'); 
  case 1
     F.forward = 1;
     F.inverse = 0;
     mol = varargin{1}; 
     if ( isa(mol,'Molecule') )
     %
     % the input arguement is a Molecule object
     %
        n1 = get(mol,'n1');
        n2 = get(mol,'n2');
        n3 = get(mol,'n3');
        F.n1 = n1;
        F.n2 = n2;
        F.n3 = n3;
        grid = Ggrid(mol);
        F.vol = det(get(mol,'supercell'));
        F.idxnz = get(grid,'idxnz');
        F = class(F,'KSFFT'); 
     else 
        error('The input must be a Molecule object');
     end;
  otherwise
     error('Must have one arguement');
end;
