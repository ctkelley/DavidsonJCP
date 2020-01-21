clear all;
diary runC60_davidson.out
ecuts = [6.25 12.5 25];
%ecuts = [6.25];
opts  = setksopt;
opts.maxscfiter = 1;
opts.maxcgiter = 50;
opts.eigmethod = 'davidson';
opts.verbose = 'on';
for j = 1:length(ecuts)
   fprintf('\n');
   fprintf('Running Ecut = %11.3e Hartree\n', ecuts(j));
   C60_setup;
   mol = set(mol,'ecut',ecuts(j));
   [mol,H,X,info]=scf(mol,opts);
   clear mol;
end;
diary off
