function display(grid)
%
% display method for atom
%
fprintf('Energy cutoff used: %11.3e\n', grid.ecut);
ng = grid.ng;
if ( iscell(ng) )
   numkpts = length(ng);
   for kpt = 1:numkpts
      fprintf('kpoint %d:\n', kpt);
      fprintf('   number of nonzeros within ecut: %d\n', grid.ng{kpt});
   end
else
   %
   % no a solid, no k point
   %
   fprintf('number of nonzeros within ecut: %d\n', grid.ng);
   fprintf('          gx             gy              gz             |g|^2\n');
   for i = 1:grid.ng
      fprintf('%15.6e %15.6e %15.6e %15.6e\n', ...
              grid.gkx(i), grid.gky(i),grid.gkz(i), grid.gkk(i));
   end
end;
