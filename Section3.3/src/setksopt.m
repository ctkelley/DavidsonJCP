function options = setksopt(varargin)
%
% Purpse: Create and/or alter KSSOLV options structure
%
% Usage: 
%     options = ksoptset('param1',value1,'param2',value2,...);
%     creates a KSSOLV option structure options in which the named 
%     parameters have the specified values.  Any unspecified parameters 
%     are set to defaults values.
%
%     newoptions = ksoptset(oldoptions,'param1',value1,'param2',value2,...);
%     resets a KSSOLV option structure oldoptions in which the named 
%     parameters have the specified values.  Any unspecified parameters 
%     are copied from oldoptions to newoptions.
%
%     The list of options and allowed values are:
%
%    option name     allowed values         purpose   
%    -----------     --------------         -------
%
%      verbose       'on' or 'off'        display diagonostic info
%     eigmethod      'lobpcg'             Method for computing the 
%                    'eigs'               the invariant subspace associated
%                                         with the lowest (occupied) states
%                                         of a fixed KS Hamiltonian
%   maxscfiter       positive integer     maximum SCF iterations allowed
%   maxdcmiter       positive integer     maximum DCM iterations allowed
%   maxinerscf       positive integer     maximum number of inner SCF 
%                                         iterations allowed in the DCM 
%                                         algorithm
%    maxcgiter       positive integer     maximum LOBPCG iterations allowed
%  maxeigsiter       positive integer     maximum implicit restarts allowed
%                                         in the eigs function
%       scftol       positive float       convergence tolerance for SCF
%       dcmtol       positive float       convergence tolerance for DCM
%        cgtol       positive float       convergence tolerance for LOBPCG
%      eigstol       positive float       convergence tolerance for eigs
%     what2mix       'pot' (default)      The object to be mixed during
%                    'rho'                the SCF iteration
%
%      mixtype       'anderson'           charge or potential mixing types
%                    'broyden'
%                    'broyden1'
%                    'pulay'
%                    'kerker'
%                    'pulay+kerker'
%                    'simple'
%                    'off'
%
%       mixdim       positive integer     maximum number of previous potentials
%                                         to be mixed
%      betamix       positive float       The damping parameter used in 
%                                         Anderson and Broyden mixing
%        brank       positive integer     The rank of the Broyden update
%
%           X0       Wavefun object       An initial guess of the wavefunctions
%         rho0       3-D array            An initial guess of the charge
%                                         To restart an SCF or DCM calculation
%                                         using previously computed results
%                                         One should pass in both the previously
%                                         obtained wavefunction and charge.
%
%       degree       positive integer     The degree of the Chebyshev polynomial
%                                         used in chebyscf
%        force       bool                 Whether calculate force
%                               


if (nargin == 0) && (nargout == 0)
    fprintf('                verbose: [ off | on ]\n');
    fprintf('              eigmethod: [ lobpcg | eigs | omm ]\n');
    fprintf('             maxscfiter: [ positive integer ]\n');
    fprintf('             maxdcmiter: [ positive integer ]\n');
    fprintf('             maxinerscf: [ positive integer ]\n');
    fprintf('              maxcgiter: [ positive integer ]\n');
    fprintf('            maxeigsiter: [ positive integer ]\n');
    fprintf('                 scftol: [ positive float ]\n');
    fprintf('                 dcmtol: [ positive float ]\n');
    fprintf('                  cgtol: [ positive float ]\n');
    fprintf('                 eigtol: [ positive float ]\n');
    fprintf('               what2mix: [ pot | rho |]\n');
    fprintf('                mixtype: [ anderson | broyden | pulay | pulay+kerker | simple | off ]\n');
    fprintf('                 mixdim: [ positive integer ]\n'); 
    fprintf('                betamix: [ positive float ]\n');
    fprintf('                  brank: [ positive integer ]\n'); 
    fprintf('                     X0: [ Wavefun object ]\n');
    fprintf('                   rho0: [ 3-D array ]\n');
    fprintf('                 degree: [ positive integer ]\n');
    fprintf('                  force: [ bool ]\n');
    fprintf('                 factorOrbitals: [ positive float ]\n');
    
    fprintf('\n');
    return;
end
% Create a struct of all the fields with all values set to 
allfields = {'verbose';'eigmethod';'maxscfiter';'maxdcmiter';...
             'maxinerscf';'maxcgiter';'maxeigsiter';...
             'scftol';'dcmtol';'cgtol';'eigstol';'what2mix';'mixtype';...
             'mixdim';'betamix';'brank';'X0';'rho0';'degree';'force';...
             'factorOrbitals'};

% create cell array
structinput = cell(2,length(allfields));
% fields go in first row
structinput(1,:) = allfields';
% []'s go in second row
structinput(2,:) = {[]};
% turn it into correctly ordered comma separated list and call struct
options = struct(structinput{:});

% set default values
options.verbose = 'off';
options.eigmethod = 'lobpcg';
options.maxscfiter = 10;
options.maxdcmiter = 10;
options.maxcgiter  = 10;
options.maxeigsiter = 300;
options.maxinerscf = 3;
options.scftol  = 1e-8;
options.dcmtol  = 1e-8;
options.cgtol   = 1e-9;
options.eigstol = 1e-10;
options.what2mix = 'pot';
options.mixtype = 'anderson';
options.mixdim = options.maxscfiter-1;
options.betamix = 0.8;
options.brank = 1;
options.degree = 10;
options.force = 1;
options.factorOrbitals = 1;

numberargs = nargin; % we might change this value, so assign it

Names = allfields;
m = size(Names,1);
names = lower(Names);

i = 1;
% check whether the first input argument is an option structure
while i <= numberargs
    arg = varargin{i};
    if ischar(arg)                         % arg is an option name
        break;
    end
    if ~isempty(arg)                      % [] is a valid options argument
        if ~isa(arg,'struct')
            error('MATLAB:setksopt:NoParamNameOrStruct',...
                ['Expected argument %d to be a string parameter name ' ...
                'or an options structure\ncreated with SETKSOPT.'], i);
        end
        for j = 1:m
            if any(strcmp(fieldnames(arg),Names{j,:}))
                val = arg.(Names{j,:});
            else
                val = [];
            end
            if ~isempty(val)
                if ischar(val)
                    val = lower(deblank(val));
                end
                checkfield(Names{j,:},val);
                options.(Names{j,:}) = val;
            end
        end
    end
    i = i + 1;
end

% A finite state machine to parse name-value pairs.
if rem(numberargs-i+1,2) ~= 0
    error('MATLAB:setksopt:ArgNameValueMismatch',...
        'Arguments must occur in name-value pairs.');
end
expectval = 0;                          % start expecting a name, not a value
while i <= numberargs
    arg = varargin{i};

    if ~expectval
        if ~ischar(arg)
            error('MATLAB:setksopt:InvalidParamName',...
                'Expected argument %d to be a string parameter name.', i);
        end

        lowArg = lower(arg);
        j = strmatch(lowArg,names);
        if isempty(j)             % if no matches
            [wasinmatlab, optionname] = checkkssolvonlylist(lowArg);
            if ~wasinmatlab
                error('MATLAB:setksopt:InvalidParamName',...
                    'Unrecognized parameter name ''%s''.', arg);
            else
                warning('MATLAB:setksopt:InvalidParamName',...
                    ['The option ''%s'' is an Optimization Toolbox option and is not\n', ...
                     'used by any MATLAB functions. This option will be ignored and not included\n', ...
                     'in the options returned by SETKSOPT. Please change your code to not use \n', ...
                     'this option as it will error in a future release.'], ...
                     optionname);
                i = i + 2; % skip this parameter and its value; go to next parameter
                continue; % skip the rest of this loop
            end
        elseif length(j) > 1                % if more than one match
            % Check for any exact matches (in case any names are subsets of others)
            k = strmatch(lowArg,names,'exact');
            if length(k) == 1
                j = k;
            else
                msg = sprintf('Ambiguous parameter name ''%s'' ', arg);
                msg = [msg '(' Names{j(1),:}];
                for k = j(2:length(j))'
                    msg = [msg ', ' Names{k,:}];
                end
                msg = sprintf('%s).', msg);
                error('MATLAB:setksopt:AmbiguousParamName', msg);
            end
        end
        expectval = 1;                      % we expect a value next

    else
        if ischar(arg)
            arg = lower(deblank(arg));
        end
        checkfield(Names{j,:},arg);
        options.(Names{j,:}) = arg;
        expectval = 0;
    end
    i = i + 1;
end

function [wasinmatlab, optionname] = checkkssolvonlylist(lowArg)
% Check if the user is trying to set an option that is only used by
% KSSOLV functions -- this used to have no effect.
% Now it will warn. In a future release, it will error.  
names =  {...
    'verbose'; ...
    'eigmethod'; ...
    'maxscfiter'; ...
    'maxdcmiter'; ...
    'maxinerscf'; ...
    'maxcgiter'; ...
    'maxeigsiter'; ...
    'scftol'; ...
    'dcmtol'; ...
    'cgtol'; ...
    'eigstol'; ...
    'what2mix'; ...
    'mixtype'; ...
    'mixdim'; ...
    'betamix'; ...
    'brank'; ...
    'X0'; ...
    'rho0'; ...
    'degree';
    'force'};
lowernames = lower(names);
k = strmatch(lowArg,lowernames);
wasinmatlab = ~isempty(k);
if wasinmatlab
    optionname = names{k};
else
    optionname = '';
end

%--------------------------------------------------------------
function checkfield(field,value)
%CHECKFIELD Check validity of structure field contents.
%   CHECKFIELD('field',V) checks the contents of the specified
%   value V to be valid for the field 'field'. 
%

% empty matrix is always valid
if isempty(value)
    return
end

validfield = true;
switch field
    case {'verbose'} % off,on
        [validvalue, errmsg, errid] = onOffType(field,value);
    case {'eigmethod'} % lobpcg,eigs,omm
        [validvalue, errmsg, errid] = eigType(field,value);
    case {'scftol','dcmtol','cgtol','eigstol','betamix','factorOrbitals'} % real scalar
        [validvalue, errmsg, errid] = nonNegReal(field,value);
    case {'maxscfiter','maxdcmiter','maxinerscf','maxcgiter','maxeigsiter','brank','mixdim'} % integer 
        [validvalue, errmsg, errid] = nonNegInteger(field,value);
    case {'what2mix'} % string for mixing types
        [validvalue, errmsg, errid] = mixObj(field,value);
    case {'mixtype'} % string for mixing types
        [validvalue, errmsg, errid] = mixType(field,value);
    case {'X0'} % Wavefun object
        [validvalue, errmsg, errid] = WaveType(field,value);
    case {'rho0'}  %string for mixing types
        [validvalue, errmsg, errid] = ArrayType(field,value);
    case {'degree'}  %integer
        [validvalue, errmsg, errid] = nonNegInteger(field,value);
    case {'force'}  %boolean integer
        [validvalue, errmsg, errid] = nonNegInteger(field,value);
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

%------------------------------------------------------------------

function [valid, errmsg, errid] = nonNegReal(field,value)
% Any nonnegative real scalar or sometimes a special string
valid =  isreal(value) && isscalar(value) && (value >= 0) ;

if ~valid
    if ischar(value)
        errid = 'MATLAB:funfun:setksopt:NonNegReal:negativeNum';
        errmsg = sprintf('Invalid value for OPTIONS parameter %s: must be a real non-negative scalar (not a string).',field);
    elseif (value < 0)
        errid = 'MATLAB:funfun:setksopt:NonNegReal:negativeNum';
        errmsg = sprintf('Invalid value for OPTIONS parameter %s: must be a real non-negative scalar.',field);
    end
else
    errid = '';
    errmsg = '';
end

%--------------------------------------------------------------------

function [valid, errmsg, errid] = nonNegInteger(field,value,string)
% Any nonnegative real integer scalar or sometimes a special string
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

%-----------------------------------------------------------------------------

function [valid, errmsg, errid] = onOffType(field,value)
% One of these strings: on, off
valid =  ischar(value) && any(strcmp(value,{'on';'off'}));
if ~valid
    errid = 'MATLAB:funfun:setksopt:onOffType:notOnOffType';
    errmsg = sprintf('Invalid value for OPTIONS parameter %s: must be ''off'' or ''on''.',field);
else
    errid = '';
    errmsg = '';
end

%------------------------------------------------------------------------
function [valid, errmsg, errid] = mixType(field,value)
% check for mixing type
valid =  ischar(value);
if ~valid
    errid = 'MATLAB:funfun:setksopt:mixType';
    errmsg = sprintf('Invalid value for OPTIONS parameter %s: must be one of the following strings: anderson,broyden,broyden1,pulay,kerker,pulay+kerker','simple','off',field);
else
    if (any(strcmp(value,...
       {'anderson';'broyden';'broyden1';'pulay';'kerker';'pulay+kerker';...
        'Anderson';'Broyden';'Broyden1';'Pulay';'Kerker';'Pulay+Kerker';
        'simple';'Simple';'off';'Off'})))
       errid = '';
       errmsg = '';
    else
       valid = 0;
       errid = 'MATLAB:funfun:setksopt:mixType';
       errmsg = sprintf('Invalid value for OPTIONS parameter %s: must be one of the following: anderson, broyden, broyden1, pulay, kerker, pulay+kerker, simple, off\n',field);

    end
end
%------------------------------------------------------------------------
function [valid, errmsg, errid] = mixObj(field,value)
% check for mixing type
valid =  ischar(value);
if ~valid
    errid = 'MATLAB:funfun:setksopt:mixObj';
    errmsg = sprintf('Invalid value for OPTIONS parameter %s: must be one of the following strings: pot,rho',field);
else
    if (any(strcmp(value,{'pot';'rho'})))
       errid = '';
       errmsg = '';
    else
       valid = 0;
       errid = 'MATLAB:funfun:setksopt:mixObj';
       errmsg = sprintf('Invalid value for OPTIONS parameter %s: must be one of the following: pot, rho\n',field);
    end
end
%------------------------------------------------------------------------
function [valid, errmsg, errid] = eigType(field,value)
% check for mixing type
valid =  ischar(value);
if ~valid
    errid = 'MATLAB:funfun:setksopt:eigType';
    errmsg = sprintf('Invalid value for OPTIONS parameter %s: must be one of the following strings: lobpcg,eigs',field);
else
    if (any(strcmp(lower(value),{'lobpcg';'eigs';'chebyfilt';'omm'})))
       errid = '';
       errmsg = '';
    else
       valid = 0;
       errid = 'MATLAB:funfun:setksopt:eigType';
       errmsg = sprintf('Invalid value for OPTIONS parameter %s: must be one of the following: lobpcg, eigs, chebyfilt, omm\n',field);
    end
end
%-------------------------------------------------------------------------------
function [valid, errmsg, errid] = WaveType(field,value)
% check for mixing type
valid =  isa(value,'Wavefun');
if ~valid
    errid = 'MATLAB:funfun:setksopt:WaveType';
    errmsg = sprintf('Invalid value for OPTIONS parameter %s: must be a Wavefun object',field);
else
    errid = '';
    errmsg = '';
end
%------------------------------------------------------------------------
function [valid, errmsg, errid] = ArrayType(field,value)
% check for mixing type
valid =  isnumeric(value);
if ~valid
    errid = 'MATLAB:funfun:setksopt:ArrayType';
    errmsg = sprintf('Invalid value for OPTIONS parameter %s: must be a 3-D array',field);
else
    errid = '';
    errmsg = '';
end
