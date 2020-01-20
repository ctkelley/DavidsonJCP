function options = setmdopt(mol,ksopts,varargin)
% SETMDOPT Set option structure for molecule dynamics
%     options = SETMDOPT(mol,ksopts,'name1','value1','name2','value2'...)
%     creates molucule dynamics options by giving structure field name and value
%     pair or set automatically with default vaule. The structure returned by 
%     this function includes imformation of optimization method in minimizing 
%     total energy in molecule by moving atoms and some parameters required in
%     using such optimization method.
%
%     The list of options and allowed values are:
%     
%     option name     allowed values         purpose   
%     -----------     --------------         -------
%       method          1,2,3,4,5         indicating optimization method
%                                         in molecule dynamics, 1 represents
%                                         using matlab toolbox fminunc, 2 
%                                         represents using nonlinear CG by Amartya
%                                         3 as nonlinear CG by M. Overton, 4
%                                         as  BFGS, 5 as quench-by-fire method
%                                         
%       md_iter       positive integer    indicate iteration times used in
%                                         melecule, can be omitted, also 
%                                         assigned in each method
%                                         
%      int_cord          logical          index of using internal coordinate
%                                         of atoms in a molecule
%                                        
%      fminunc_opts     structure         used by function 'fminunc', consisting
%                                         field of Algorithm Diagnostics FiniteDifferenceStepSize
%                                         SpecifyObjectiveGradient MaxFunctionEvaluations
%                                         MaxIterations OptimalityTolerance
%                                         see also fminunc
%                        
%       NLCG_opts       structure         required in function 'NLCG' indicating
%                                         parameters of NLCG iteration
%                                         see also NLCG
%                                                
%       nlcg_pars       structure         required in function 'nlcg' indicating
%       nlcg_opts                         parameters and options of NLCG iteration
%                                         see also nlcg
%                                         
%       BFGS_pars       structure         required in function 'bfgs' indicating
%       BFGS_opts                         parameters and options of bfgs iteration
%                                        see also bfgs
%                                         
%       fire_pars       structure         required in function 'quenchbyfire' 
%       fire_opts                         indicating parameters and options 
%                                         of quench-fire iteration,quenchbyfire
%
% Copyright (c) 2018-2019 Bichen Lu and Yingzhou Li,
%                         Fudan University and Duke University

% Create a struct of all the fields with all values set to 
allfields={'method';'md_iter';'int_cord';'fminunc_opts';'NLCG_opts';...
'nlcg_pars';'nlcg_opts';'BFGS_pars';'BFGS_opts';'fire_pars';'fire_opts'};
structinput = cell(2,length(allfields));
structinput(1,:) = allfields';
structinput(2,:) = {[]};
options = struct(structinput{:});
options.method = 5;
options.md_iter=100;
options.int_cord= true;

xyzlist = mol.xyzlist;
x0 = xyzlist(:);
natoms=sum(mol.natoms);


Names=allfields;
i=1;
%set values to fields if values are given
while i+2<= nargin
    arg=varargin{i};
    if ischar(arg)
        break;
    end
    if ~isempty(arg)                      % [] is a valid options argument
        if ~isa(arg,'struct')
            error('MATLAB:mdopt:NoParamNameOrStruct',...
                ['Expected argument %d to be a string parameter name ' ...
                'or an options structure\ncreated with MDOPT.'], i);
        end
        for j = 1:m
            if any(strcmp(fieldnames(arg),Names{j,:}))
                val = arg.(Names{j,:});
            else
                val = [];
            end
            if ~isempty(val)
                checkfield(Names{j,:},val);
                options.(Names{j,:}) = val;
            end
        end
    end
    i=i+1;
end
if rem(nargin-i+1,2) ~= 0
    error('MATLAB:mdopt:ArgNameValueMismatch',...
        'Arguments must occur in name-value pairs.');
end
expectval = 0;   

%match input field name and field in setmdopt.m
while  i+2<= nargin
    arg = varargin{i};
    if ~expectval
        if ~ischar(arg)
            error('MATLAB:mdopt:InvalidParamName',...
                'Expected argument %d to be a string parameter name.', i);
        end

        j = strmatch(arg,Names);
        if isempty(j)             % if no matches
            error('MATLAB:setksopt:InvalidParamName',...
                'Unrecognized parameter name ''%s''.', arg);
            i = i + 2; % skip this parameter and its value; go to next parameter
            continue; % skip the rest of this loop
        elseif length(j) > 1                % if more than one match
            % Check for any exact matches (in case any names are subsets of others)
            k = strmatch(arg,Names,'exact');
            if length(k) == 1
                j = k;
            else
                msg = sprintf('Ambiguous parameter name ''%s'' ', arg);
                msg = [msg '(' Names{j(1),:}];
                for k = j(2:length(j))'
                    msg = [msg ', ' Names{k,:}];
                end
                msg = sprintf('%s).', msg);
                error('MATLAB:mdopt:AmbiguousParamName', msg);
            end
        end
        expectval = 1;                      % we expect a value next

    else
        if ischar(arg)
            arg = deblank(arg);
        end
        checkfield(Names{j,:},arg);
        options.(Names{j,:}) = arg;
        expectval = 0;
    end
    i = i + 1;
end
switch options.method
    case 1
%         set default values for using toolbox fminunc, the setting will be done 
%         through optimoptions function or set manually, note different 
%         version of fminunc requires different setting. Important setting 
%         includes optimization algorithm, max iteration number, tolerance.
        if isempty(options.fminunc_opts)
            options.fminunc_opts=optimoptions('fminunc');
            options.fminunc_opts.Algorithm = 'quasi-newton';
            options.fminunc_opts.Diagnostics = 'on';
            if verLessThan('optim', '7.4')
                options.fminunc_opts.FinDiffRelStep = 0.01;
                options.fminunc_opts.GradObj = 'on';
                options.fminunc_opts.MaxFunEvals=100;
                options.fminunc_opts.MaxIter = 100;
                options.fminunc_opts.TolX = 5e-3/sqrt(3*natoms);
            else
                options.fminunc_opts.FiniteDifferenceStepSize = 0.01;
                options.fminunc_opts.SpecifyObjectiveGradient = true;
                options.fminunc_opts.MaxFunctionEvaluations = 100;
                options.fminunc_opts.MaxIterations = 100;
                options.fminunc_opts.OptimalityTolerance = 5e-3/sqrt(3*natoms);
            end
        end
    case 2
%         set default values for NLCG options, NLCG requires options of outer
%         and inner iteration times, outer and inner iteration tolerance and 
%         baisc step size setting
        if isempty(options.NLCG_opts)
            NLCG_field = {'outiter_max';'ineriter_max';'NLCG_tol1';'NLCG_tol2';'step_size0'};
            NLCG_opts = cell(2,length(NLCG_field));
            NLCG_opts(1,:) = NLCG_field';
            NLCG_opts(2,:) = {[]};
            options.NLCG_opts = struct(NLCG_opts{:});
            options.NLCG_opts.outiter_max = 100;
            options.NLCG_opts.ineriter_max = 6;
            options.NLCG_opts.n = 10;
            options.NLCG_opts.NLCG_tol1 = 1e-6;
            options.NLCG_opts.NLCG_tol2 = 1e-8;
            options.NLCG_opts.step_size0 = 0.02;
        end
    case 3
%         set default values for nlcg options, nlcg requires two inputs, pars
%         and opts, pars including nvar(number of variable dimension), fdname
%         refering to a m-file calculating function value and gradient, int_cord
%         , mol, ksopts. opts includes start point of iterration x0, maxit
%         (max times of iteration), normtol(tolrance of gradient normal).
%         strongwolfe, wolfe1 and wolfe2 refering to parmaters in wolve condition 
%         in line search, version is chosen as follow
        % options.version (used to obtain different choices of beta):
        % 'P' for Polak-Ribiere-Polyak (not recommended: fails on hard problems)
        % 'F' for Fletcher-Reeves (not recommended: often stagnates)
        % 'C' for Polak-Ribiere-Polyak Constrained by Fletcher-Reeves
        % (recommended, combines advantages of 'P' and 'F'; default)
        % 'S' for Hestenes-Stiefel (not recommended)
        % 'Y' for Dai-Yuan (allows weak Wolfe line search, see nlcg.m)
        % 'Z' for Hager-Zhang
        % '-' for Steepest Descent (for comparison)
        
        if isempty(options.nlcg_pars)
            nlcg_pars.nvar = length(mol.xyzlist(:));
            nlcg_pars.fgname = 'ksef';
            nlcg_pars.int_cord = true;
            nlcg_pars.mol = mol;
            nlcg_pars.ksopts = ksopts;
            options.nlcg_pars = nlcg_pars;
        end
        if isempty(options.nlcg_opts)
            nlcg_options.x0 = x0;
            nlcg_options.nstart = 1;
            nlcg_options.maxit = 100;
            nlcg_options.normtol = 1.0e-4;
            nlcg_options.version = 'C';
            nlcg_options.strongwolfe = 1;
            nlcg_options.wolfe1 = 0.1;           % 0 < c1 < c2 < 1/2
            nlcg_options.wolfe2 = 0.49;          % c2 < 1/2 for NLCG
            nlcg_options.prtlevel = 1;
            options.nlcg_opts = nlcg_options;
        end
    case 4
%         set default values for BFGS options, BFGS requires two inputs, pars
%         and opts, mainly same as nlcg
        if isempty(options.BFGS_pars) 
            BFGS_pars.nvar = length(x0);
            BFGS_pars.fgname = 'ksef'; 
            BFGS_pars.mol = mol;
            BFGS_pars.ksopts = ksopts;
            options.BFGS_pars=BFGS_pars;
        end
        if isempty(options.BFGS_opts) 
            BFGS_opts.x0 = x0;
            BFGS_opts.nstart = 1;
            BFGS_opts.maxit = 1000;
            BFGS_opts.nvec = 0;           % 0 => full BFGS, else specify m \in [3, 20]
            BFGS_opts.H0 = sparse(eye(BFGS_pars.nvar));
            BFGS_opts.scale = 1;
            BFGS_opts.normtol = 1.0e-6;   % default
            BFGS_opts.strongwolfe = 1;
            BFGS_opts.wolfe1 = 1.0e-4;    % 0 < c1 < c2 < 1
            BFGS_opts.wolfe2 = 0.5;       % c2 = 1/2 is also default
            BFGS_opts.ngrad = 1;          % saves the final gradient
            BFGS_opts.prtlevel = 1;       % default
            options.BFGS_opts = BFGS_opts;
        end
    case 5
%         set default values for BFGS options, BFGS requires two inputs, pars
%         and opts, mainly same as previous setting, dt refering to initial
%         setting of step size, recommanded 1 femtosecond, equal to 41.34
%         in atomic time unit
        if isempty(options.fire_pars) 
            fire_pars.mass = 4.0;
            fire_pars.fDec = 0.5;
            fire_pars.fInc = 1.1;
            fire_pars.nMin = 5;
            fire_pars.alphaStart = 0.1;
            fire_pars.fAlpha = 0.99;
            options.fire_pars = fire_pars;
        end
        if isempty(options.fire_opts)
            fire_opts.dt = 41.3413745758;
            fire_opts.MAXITER = 100;
            fire_opts.TOL = 1.0E-4;
            options.fire_opts=fire_opts;
        end
end





function checkfield(field,value)

% empty matrix is always valid
if isempty(value)
    return
end

validfield = true;
switch field
    case {'md_iter'} % real scalar
        [validvalue, errmsg, errid] = nonNegInteger(field,value);
    case {'method'} % integer 
        [validvalue, errmsg, errid] = mdmethod(field,value);
    case {'int_cord'} % string for mixing types
        [validvalue, errmsg, errid] = logic(field,value);
    case {'fminunc_opts','NLCG_opts',...
    'nlcg_pars','nlcg_opts','BFGS_pars','BFGS_opts','fire_pars','fire_opts'} % string for mixing types
        [validvalue, errmsg, errid] = methodstruct(field,value);
    otherwise
        validfield = false;  
        validvalue = false;
        errmsg = sprintf('Unrecognized parameter name ''%s''.', field);
        errid = 'MATLAB:setksopt:checkfield:InvalidParamName';
end

if validvalue 
    return;
else
    % Throw the MATLAB invalid value error
    error(errid, errmsg);
end


function [valid, errmsg, errid] = logic(field,value,string)
% Any nonnegative real integer scalar or sometimes a special string
valid =  islogical(value) ||(value==1)||(value==0);
if ~valid
    errid = 'MATLAB:funfun:setksopt:logic:not a logical var';
    errmsg = sprintf('Invalid value for OPTIONS parameter %s: must be a logical var or 0/1 var.',field);

else
    errid = '';
    errmsg = '';
end


function [valid, errmsg, errid] = nonNegInteger(field,value,string)

valid =  isreal(value) && isscalar(value) && (value >= 0) && value == floor(value) ;
if ~valid
    if ischar(value)
        errid = 'MATLAB:funfun:setksopt:nonNegInteger:notANonNegInteger';
        errmsg = sprintf('Invalid value for OPTIONS parameter %s: must be a real non-negative integer (not a string).',field);
    elseif (value < 0) 
        errid = 'MATLAB:funfun:setksopt:nonNegInteger:notANonNegInteger';
        errmsg = sprintf('Invalid value for OPTIONS parameter %s: must be a real non-negative integer.',field);
    end
else
    errid = '';
    errmsg = '';
end


function [valid, errmsg, errid] = mdmethod(field,value,string)

valid =  (value==1)||(value==2)||(value==3)||(value==4)||(value==5) ;
if ~valid
    errid = 'MATLAB:funfun:setksopt:mdmethod:notAmdmethod';
    errmsg = sprintf('Invalid value for OPTIONS parameter %s: must be a choice among {1,2,3,4,5} to represent a optimization method.',field);
else
    errid = '';
    errmsg = '';
end


function [valid, errmsg, errid] = methodstruct(field,value,string)
valid =  isstruct(value) ;
if ~valid
    errid = 'MATLAB:funfun:setksopt:methodstruct:notAStruct';
    errmsg = sprintf('Invalid value for OPTIONS parameter %s: must be a struct to set the paramaters of the md_method.',field);
else
    errid = '';
    errmsg = '';
end
