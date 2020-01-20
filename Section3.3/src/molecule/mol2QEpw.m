function filename = mol2QEpw(mol,filename,method,ksopt)
% MOL2QEPW generate input file for Quantum-Espesso pw.x.
%   MOL2QEPW(mol,filename,method,ksopt) generate the input file for
%   Quantum-Espesso pw.x based on molecule mol and the given method and
%   KSSOLV option ksopt. And the input file is saved as filename. By
%   default, the filename is [mol.name '.in'], the method is scf, and the
%   ksopt is the default ksopt.
%
%   See also Molecule, Crystal, setksopt.

%  Copyright (c) 2016-2017 Yingzhou Li and Chao Yang,
%                          Stanford University and Lawrence Berkeley
%                          National Laboratory
%  This file is distributed under the terms of the MIT License.

if nargin < 4 || isempty(ksopt)
    ksopt = setksopt();
end
if nargin < 3 || isempty(method)
    method = 'scf';
end
molname = strrep(mol.name,' ','_');
if nargin < 2 || isempty(filename)
    filename = [molname '.in'];
end

fp = fopen(filename,'w');

%=========================================================================
%                                 CONTROL
%-------------------------------------------------------------------------
% Setting for Control
if strcmpi(ksopt.verbose,'on')
    verbose = 'high';
else
    verbose = 'low';
end
if ksopt.force
    force = '.TRUE.';
else
    force = '.FALSE.';
end
ppdir = [ kssolvroot 'ppdata/' ];
if isempty(kssolvpptype) || strcmpi(kssolvpptype,'default')
    ppdir = [kssolvroot 'ppdata/default/'];
end
%-------------------------------------------------------------------------
% print CONTROL part
fprintf(fp,'&CONTROL\n');
fprintf(fp,'    calculation  = ''%s'',\n',method);
fprintf(fp,'    title        = ''%s'',\n',mol.name);
fprintf(fp,'    verbosity    = ''%s'',\n',verbose);
fprintf(fp,'    restart_mode = ''from_scratch'',\n');
fprintf(fp,'    tprnfor      = %s,\n',force);
fprintf(fp,'    outdir       = ''./tmp'',\n');
fprintf(fp,'    prefix       = ''%s'',\n',molname);
fprintf(fp,'    pseudo_dir   = ''%s'',\n',ppdir);
fprintf(fp,'/\n\n');
%=========================================================================


%=========================================================================
%                                 SYSTEM
%-------------------------------------------------------------------------
% Setting for System
if mol.temperature > 0
    occupations = 'smearing';
    % Kelvin to Rydberg
    degauss = mol.temperature*6.336857346553284e-06;
    smearing = 'fermi-dirac';
else
    occupations = 'fixed';
    degauss = 0;
    smearing = 'gaussian';
end

%-------------------------------------------------------------------------
% print SYSTEM part
fprintf(fp,'&SYSTEM\n');
fprintf(fp,'    ibrav        = 0,\n');
fprintf(fp,'    nat          = %d,\n',length(mol.alist));
fprintf(fp,'    ntyp         = %d,\n',length(mol.atoms));
fprintf(fp,'    ecutwfc      = %f,\n',mol.ecut*2*meDef());
fprintf(fp,'    occupations  = ''%s'',\n',occupations);
fprintf(fp,'    degauss      = %.10f,\n',degauss);
fprintf(fp,'    smearing     = ''%s'',\n',smearing);
fprintf(fp,'    nspin        = %d,\n',mol.nspin);
% TODO: KSSOLV only support PZ-LDA
% fprintf(fp,'    input_dft    = ''pz''\n');
fprintf(fp,'/\n\n');
%=========================================================================


%=========================================================================
%                                ELECTRONS
%-------------------------------------------------------------------------
% Setting for Electrons
if strcmpi(method,'scf')
    scfmaxiter = ksopt.maxscfiter;
    convthr = ksopt.scftol;
else
    scfmaxiter = ksopt.maxinerscf;
    convthr = ksopt.dcmtol;
end
if strcmpi(ksopt.eigmethod,'lobpcg')
    diagmethod = 'cg';
    diag_thr = ksopt.cgtol;
    cgmaxiter = ksopt.maxcgiter;
elseif strcmpi(ksopt.eigmethod,'eig')
    diagmethod = 'david';
    diag_thr = ksopt.eigstol;
    cgmaxiter = ksopt.maxcgiter;
else
    % TODO: find OMM option in QE
    diagmethod = 'david';
    diag_thr = ksopt.eigstol;
    cgmaxiter = ksopt.maxcgiter;
end
%-------------------------------------------------------------------------
% print ELECTRONS part
fprintf(fp,'&ELECTRONS\n');
fprintf(fp,'    electron_maxstep = %d,\n',scfmaxiter);
fprintf(fp,'    conv_thr         = %d,\n',convthr);
% TODO: Match the mixing mode with QE mixing
fprintf(fp,'    mixing_mode      = ''plain'',\n');
fprintf(fp,'    mixing_beta      = %f,\n',ksopt.betamix);
fprintf(fp,'    mixing_ndim      = %d,\n',ksopt.mixdim);
fprintf(fp,'    diagonalization  = ''%s'',\n',diagmethod);
fprintf(fp,'    diago_thr_init   = %e,\n',diag_thr);
fprintf(fp,'    diago_cg_maxiter = %d,\n',cgmaxiter);
fprintf(fp,'/\n\n');
%=========================================================================


%=========================================================================
%                              ATOMIC_SPECIES
%-------------------------------------------------------------------------
% print ATOMIC_SPECIES
fprintf(fp,'ATOMIC_SPECIES\n');
for it = 1:length(mol.atoms)
    asym = mol.atoms(it).symbol;
    amass = mol.atoms(it).amass;
    
    [pptype,ppext] = kssolvpptype();
    
    ppregexp = [asym '[._-]' pptype '[\w.-]*(?i)(' ppext ')(?-i)'];
    pppath = [kssolvroot 'ppdata/'];
    [~,ppfile] = dirregexp(pppath,ppregexp);
    if isempty(ppfile) || strcmpi(pptype,'default')
        pppath = [kssolvroot 'ppdata/default/'];
        ppregexp = [asym '[._-][\w.-]*(?i)(' ppext ')(?-i)'];
        [~,ppfile] = dirregexp(pppath,ppregexp);
    end
    
    if isempty(ppfile)
        warning('The pseudopotential files are not specified.');
    end
    fprintf(fp, '%3s %18.10e   %s\n', asym, amass, ppfile);
end
fprintf(fp,'\n');
%=========================================================================


%=========================================================================
%                             ATOMIC_POSITIONS
%-------------------------------------------------------------------------
% print ATOMIC_POSITIONS
fprintf(fp,'ATOMIC_POSITIONS {bohr}\n');
for it = 1:length(mol.alist)
    asym = mol.atoms(mol.alist(it)).symbol;
    xyz = mol.xyzlist(it,:);
    fprintf(fp, '%3s %18.10e %18.10e %18.10e\n', asym, xyz);
end
fprintf(fp,'\n');
%=========================================================================


%=========================================================================
%                                K_POINTS
%-------------------------------------------------------------------------
% print K_POINTS
if isa(mol,'Crystal')
    fprintf(fp,'K_POINTS tpiba\n');
    fprintf(fp,'%4d\n',mol.nkpts);
    kpts = mol.kpts/2/pi*mol.supercell;
    for it = 1:mol.nkpts
        fprintf(fp, '%18.10e %18.10e %18.10e %18.10e\n', ...
            kpts(it,:), mol.wks(it));
    end
    fprintf(fp,'\n');
end
%=========================================================================


%=========================================================================
%                               CELL_PARAMETERS
%-------------------------------------------------------------------------
% print CELL_PARAMETERS
fprintf(fp,'CELL_PARAMETERS {bohr}\n');
for it = 1:3
    fprintf(fp,' %18.10e   %18.10e   %18.10e\n', mol.supercell(it,:) );
end
fprintf(fp,'\n');
%=========================================================================

fclose(fp);
end