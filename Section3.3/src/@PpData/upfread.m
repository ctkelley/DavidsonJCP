function pp = upfread( pp, filename )
% UPFREAD  Parse an unified pseudopotential format (UPF) document into
%          PpData type.
%    pp = UPFREAD(pp,filename) reads a file in the string input argument
%    filename.  The function assign the values in the file to the pp which
%    is of type PpData.  If the version is older than 2.0.1, the pp will be
%    read from upfreadold function.
%
% 
%    See also xmlread, PpData, upfreadold.
%
%    Reference page in Help browser
%       doc upfread

%  Copyright (c) 2016-2017 Yingzhou Li and Chao Yang,
%                          Stanford University and Lawrence Berkeley
%                          National Laboratory
%  This file is distributed under the terms of the MIT License.

e2 = e2Def();

%-------------------------------------------------------------------------
% Read the file into a long string.
str = fileread(filename);

%-------------------------------------------------------------------------
% Remove all comments in the upf file.
str = regexprep( str, '<\-.*?\->', '' ) ;


%-------------------------------------------------------------------------
% Read the version of upf file.
% In the case the version is not specified, empty string is assigned.
[UPFstr,UPFattr] = upfreadfield(str,'UPF');

pp.info.pptype = 'upf';
pp.info.version = mapvalue(UPFattr{1},'version','');
if ~strcmp(pp.info.version,'2.0.1')
    warning('upf version is not 2.0.1.');
    pp = upfreadold( pp, filename);
    return;
end
UPFstr = cell2val(UPFstr,str);


%-------------------------------------------------------------------------
% Read the PP_INFO of upf file.
% If there are multiple PP_INFO fields exist, they will be merged into a
% single string and stored in upf.info.
[infostr,~] = upfreadfield(UPFstr,'PP_INFO');

pp.info.message = '';
for infos = infostr
    pp.info.message = strcat(pp.info.message,cell2val(infos));
end


%-------------------------------------------------------------------------
% Read the PP_HEADER of upf file.
% If there are multiple PP_HEADER fields exist, it is invalid.
headerkey = {'generated','author','date','comment','element',...
             'pseudo_type','relativistic','functional','is_ultrasoft',...
             'is_paw','is_coulomb','has_so','has_wfc','has_gipaw',...
             'paw_as_gipaw','core_correction','z_valence',...
             'total_psenergy','wfc_cutoff','rho_cutoff','l_max',...
             'l_max_rho','l_local','mesh_size','number_of_wfc',...
             'number_of_proj'};
headerinit = {'','','','','',...
             '','','',false,...
             false,false,false,false,false,...
             false,false,NaN,...
             NaN,NaN,NaN,NaN,...
             NaN,NaN,NaN,NaN,...
             NaN};

[~,attr] = upfreadfield(UPFstr,'PP_HEADER');

pp.info = upfassignattr(pp.info,attr{1},headerkey,headerinit);

upf_mesh_size = pp.info.mesh_size;
upf_nbeta = pp.info.number_of_proj;

pp.anum = sym2num(pp.info.element);
pp.venum = pp.info.z_valence;

%-------------------------------------------------------------------------
% Read the PP_MESH of upf file.
% If there are multiple PP_MESH fields exist, it is invalid.
meshkey = {'dx','mesh','xmin','rmax','zmesh'};
meshinit = {NaN,NaN,NaN,NaN,NaN};

[meshstr,attr] = upfreadfield(UPFstr,'PP_MESH');
pp.info = upfassignattr(pp.info,attr{1},meshkey,meshinit);

% Read PP_R from PP_MESH
[Rstr,~] = upfreadfield(meshstr,'PP_R');
pp.r = upfreadtable(Rstr);

% Read PP_RAB from PP_MESH
[RABstr,~] = upfreadfield(meshstr,'PP_RAB');
pp.rab = upfreadtable(RABstr);


%-------------------------------------------------------------------------
% Read the PP_NLCC of upf file.
% If there are multiple PP_NLCC fields exist, it is invalid. However, the
% field itself is optional.
[nlccstr,~] = upfreadfield(UPFstr,'PP_NLCC');

pp.rho_atc = upfreadtable(nlccstr);


%-------------------------------------------------------------------------
% Read the PP_SEMILOCAL of upf file.
% If there are multiple PP_SEMILOCAL fields exist, it is invalid. However,
% the field itself is optional. The original input is in the unit of
% Ry^{1/2}.
semilocalkey = {'L','J'};
semilocalinit = {NaN,NaN};

[semilocalstr,~] = upfreadfield(UPFstr,'PP_SEMILOCAL');

semilocalstr = cell2val(semilocalstr);
if ~isempty(semilocalstr)
    pp.semilocal = struct;
    pp.semilocal = upfinitattrtable(pp.semilocal,semilocalkey,...
                                     semilocalinit,upf_nbeta);
    pp.semilocal.VNL = zeros(upf_mesh_size,upf_nbeta);
    
    semilocalstr = cell2val(semilocalstr,'');
    for ibeta = 1:upf_nbeta
        [vnlstr,attr] = upfreadfield(semilocalstr,...
                                     ['PP_VNL' num2str(ibeta)]);

        pp.semilocal = upfassignattr(pp.semilocal,attr{1},semilocalkey,...
                                     semilocalinit,ibeta);
        pp.semilocal.VNL(:,ibeta) = upfreadtable(vnlstr);
    end
    
    % convert to current unit
    pp.semilocal.VNL = pp.semilocal.VNL*sqrt(e2/2);
end


%-------------------------------------------------------------------------
% Read the PP_LOCAL of upf file.
% If there are multiple PP_LOCAL fields exist, it is invalid. The original
% input is in the unit of Ry.
[localstr,~] = upfreadfield(UPFstr,'PP_LOCAL');

pp.vloc = upfreadtable(localstr);

% convert to current unit
pp.vloc = pp.vloc/2*e2;


%-------------------------------------------------------------------------
% Read the PP_NONLOCAL of upf file.
% If there are multiple PP_NONLOCAL fields exist, it is invalid. The
% original input, D is in the unit of Ry; and beta is in the unit of
% Bohr^{1/2}.
nonlocbetakey = {'label','angular_momentum','cutoff_radius_index',...
                 'cutoff_radius','ultrasoft_cutoff_radius',...
                 'norm_conserving_radius'};
nonlocbetainit = {'',NaN,NaN,...
                  NaN,NaN,...
                  NaN};

[nonlocstr,~] = upfreadfield(UPFstr,'PP_NONLOCAL');

nonlocstr = cell2val(nonlocstr);
if ~isempty(nonlocstr)
    pp.nonloc = struct;
    pp.nonloc = upfinitattrtable(pp.nonloc,nonlocbetakey,...
                                 nonlocbetainit,upf_nbeta);
    pp.nonloc.nbeta = upf_nbeta;
    pp.nonloc.beta = zeros(upf_mesh_size,upf_nbeta);
    pp.nonloc.lll  = zeros(1,upf_nbeta);
    
    for ibeta = 1:upf_nbeta
        [betastr,attr] = upfreadfield(nonlocstr,...
            ['PP_BETA.' num2str(ibeta)]);
    
        pp.nonloc = upfassignattr(pp.nonloc,attr{1},nonlocbetakey,...
                                  nonlocbetainit,ibeta);
        pp.nonloc.beta(:,ibeta) = upfreadtable(betastr);
    end
    pp.nonloc.lll = pp.nonloc.angular_momentum;
    
    [dijstr,~] = upfreadfield(nonlocstr,'PP_DIJ');
    pp.nonloc.dij = upfreadtable(dijstr,upf_nbeta,upf_nbeta);
    
    % convert to current unit
    pp.nonloc.dij = pp.nonloc.dij/2*e2;
    
    if pp.info.is_ultrasoft || pp.info.is_paw
        % TODO
    end
end


%-------------------------------------------------------------------------
% Read the PP_RHOATOM of upf file.
% If there are multiple PP_RHOATOM fields exist, it is invalid.
[rhoatomstr,~] = upfreadfield(UPFstr,'PP_RHOATOM');

pp.rhoatom = upfreadtable(rhoatomstr);


%************************************************************************%
%                                                                        %
%                       Sub Functions in upfread                         %
%                                                                        %
%************************************************************************%

%=========================================================================
% Subfunction cell2val
% Convert the first element of any recursive cell to value.
    function val = cell2val(c,val)
        if nargin < 2
            val = {};
        end
        if isempty(c)
            return;
        end
        if iscell(c)
            val = cell2val(c{1},val);
        else
            if ~isempty(c)
                val = c;
            end
        end
    end % end of subfunc cell2val
%TODO: move this subfunction to the outside as a tool

%=========================================================================
% Subfunction str2type
% Convert the first string to the type of the second input. Currently, the
% type of the second input should be bool, number or string.
    function val = str2type(s,t)
        s = cell2val(s,'');
        if isnumeric(t)
            val = str2double(s);
        elseif islogical(t)
            val = (s == 'T') | (s == 't') | strcmp(s,'True') | ...
                  strcmp(s,'true') | (s == '1');
        elseif ischar(t)
            val = s;
        else
            error('Invalid input type for str2type');
        end
    end % end of subfunc str2type
%TODO: move this subfunction to the outside as a tool

%=========================================================================
% Subfunction mapvalue
% Find the value with the given key in the map.
    function val = mapvalue(m,k,dv)
        k = cell2val(k);
        len = size(m,1);
        val = dv;
        for it = 1:len
            if strcmpi(k,m{it,1})
                val = str2type(m{it,2},dv);
            end
        end
    end % end of subfunc mapvalue
%TODO: move this subfunction to the outside as a tool

%=========================================================================
% Subfunction upfinitattrtable
% Initialize the table for attributes
    function f = upfinitattrtable(f,keys,inits,nsize)
        if length(nsize) == 1
            nsize = repmat(nsize,size(inits));
        end
        for i = 1:length(keys)
            if ischar(inits{i})
            	f.(cell2val(keys{i})) = cell(1,nsize(i));
            else
                f.(cell2val(keys{i})) = zeros(1,nsize(i));
            end
        end
    end % end of subfunc upfinitattrtable

%=========================================================================
% Subfunction upfassignattr
% Assign attributes to upf field
    function f = upfassignattr(f,attr,keys,inits,idx)
        if nargin < 5
            for i = 1:length(keys)
                f.(cell2val(keys{i})) = mapvalue(attr,keys{i},inits{i});
            end
        else
            for i = 1:length(keys)
                if ischar(inits{i})
                    f.(cell2val(keys{i})){idx} = ...
                        mapvalue(attr,keys{i},inits{i});
                else
                    f.(cell2val(keys{i}))(idx) = ...
                        mapvalue(attr,keys{i},inits{i});
                end
            end
        end
    end % end of subfunc upfassignattr

%=========================================================================
% Subfunction upfreadtable
% Read table to upf content
    function t = upfreadtable(str,n,m)
        str = cell2val(str,'');
        if ~isempty(str)
            if nargin < 2
                t = cell2val(textscan(str,'%f'));
            elseif nargin < 3
                if n > 0
                    t = reshape(cell2val(textscan(str,'%f')),n,[]);
                else
                    t = [];
                end
            else
                if n*m > 0
                    t = reshape(cell2val(textscan(str,'%f')),n,m);
                else
                    t = NaN(n,m);
                end
            end
        else
            if nargin < 2
                t = [];
            elseif nargin < 3
                t = NaN(n,0);
            else
                t = NaN(n,m);
            end
        end
    end % end of subfunc upfreadtable

%=========================================================================
% Subfunction upfreadattr
% Read the attributions from the field head.
    function attr = upfreadattr(str)
        if isempty(str)
            attr = cell(0,2);
            return;
        end
        attrregexp = '[\w]+="[^><"]*"';
        Cstr = regexp(str,attrregexp,'match');
        lenCstr = length(Cstr);
        keyset = cell(lenCstr,1);
        valueset = cell(lenCstr,1);
        it = 1;
        for s = Cstr
            tmpstr = strsplit(s{1},'"','CollapseDelimiters',false);
            if length(tmpstr) ~= 3
                error('The attribution list is not valid.');
            end
            keyset{it} = cell2val(regexp(tmpstr{1},'[\w]+','match'));
            valueset{it} = cell2val(tmpstr{2});
            it = it+1;
        end
        attr = [keyset valueset];
    end % end of subfunc upfreadattr

%=========================================================================
% Subfunction upfreadfield
% Read the field from the file, if the field name does not exist, the
% empty str is returned in the fieldstr.
    function [fieldstr,attr] = upfreadfield(str,fieldname)
        str = cell2val(str,'');
        attrregexp = ['([\s]+[\w]+=(["][^><"]*["]|[' '''' '][^><]*[' '''' ']))*[\s]*'];
        [sIdxhead,eIdxhead] = regexp(str,['<' fieldname attrregexp '>']);
        [sIdxline,eIdxline] = regexp(str,['<' fieldname attrregexp '\/>']);
        [sIdxtail,eIdxtail] = regexp(str,['<\/' fieldname attrregexp '>']);
        
        if isempty(eIdxhead) && isempty(sIdxtail) && isempty(sIdxline)
            fieldstr{1} = '';
            attr{1} = cell(0,2);
            return;
        end
        
        sIdx = [sIdxhead;sIdxline;sIdxline;sIdxtail];
        eIdx = [eIdxhead;eIdxline;eIdxline;eIdxtail];
        vIdx = [ones(size(sIdxhead));ones(size(sIdxline))
                -ones(size(sIdxline));-ones(size(sIdxtail))];
        [~,I] = sort(sIdx);
        eIdx = eIdx(I);
        vIdx = vIdx(I);
        vIdx = cumsum(vIdx);
        if vIdx(end) ~= 0 || sum(vIdx<0)>0
            error(['Field ' fieldname ' does not paired.']);
        end
        vIdx = find(vIdx==0);
        it = length(vIdx);
        fieldstr = cell(it,1);
        attr = cell(it,1);
        for i = 1:it
            if i == 1
                ithead = 1;
            else
                ithead = vIdx(it-1)+1;
            end
            ittail = vIdx(it);
            if sIdx(ithead) == sIdx(ittail)
                fieldstr{it} = '';
                attr{it} = upfreadattr...
                    (str(sIdx(ithead)+length(fieldname)+1:eIdx(ithead)-2));
            else
                fieldstr{it} = str(eIdx(ithead)+1:sIdx(ittail)-1);
                attr{it} = upfreadattr...
                    (str(sIdx(ithead)+length(fieldname)+1:eIdx(ithead)-1));
            end
        end
        return;
        
    end % end of subfunc upfreadfield

end % end of func upfread
