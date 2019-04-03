function S = getfilelist(S,varargin)

% add suffix as part of callng function, rather than from script
if nargin>1
    %if ~isfield(S.(S.func),'load') || ~isfield(S.(S.func).load,'suffix') || isempty(S.(S.func).load.suffix{:})
    if iscell(varargin{1})
        varargin=varargin{:};
    end
    
    S.(S.func).load.suffix = varargin;
    if ~any(strcmp(S.(S.func).fname.parts,'suffix'))
        S.(S.func).fname.parts = [S.(S.func).fname.parts 'suffix'];
    end
    %else
    %    S.(S.func).load.suffix{:} = [S.(S.func).load.suffix{:} '_' varargin{:}];
    %end
end

if ~isfield(S,'func')
    S.func='T';
    S.T=S;
end

if iscell(S.path.file)
    S.path.file = S.path.file{:};
end

cd(S.path.file);
% if not specified in settings:
if ~isfield(S.(S.func),'study') || isempty(S.(S.func).study)
    S.(S.func).study= {''};
end
if ~isfield(S.(S.func).select,'groups') || isempty(S.(S.func).select.groups)
    S.(S.func).select.groups = {''};
end
if ~isfield(S.(S.func).select,'subjects') || isempty(S.(S.func).select.subjects)
    S.(S.func).select.subjects = {''};
end
if ~isfield(S.(S.func).select,'sessions') || isempty(S.(S.func).select.sessions)
    S.(S.func).select.sessions = {''};
end
if ~isfield(S.(S.func).select,'blocks') || isempty(S.(S.func).select.blocks)
    S.(S.func).select.blocks = {''};
end
if ~isfield(S.(S.func).select,'conds') || isempty(S.(S.func).select.conds)
    S.(S.func).select.conds = {''};
end
S.(S.func).subjects_in = S.(S.func).select.subjects;

% if no prefix/suffix specified
if ~isfield(S.(S.func),'load') || ~isfield(S.(S.func).load,'prefix')
    S.(S.func).load.prefix = {''};
end
if ~isfield(S.(S.func),'load') || ~isfield(S.(S.func).load,'suffix')
    S.(S.func).load.suffix = {''};
end
if ~isfield(S.(S.func),'save') || ~isfield(S.(S.func).save,'prefix')
    S.(S.func).save.prefix = {''};
end
if ~isfield(S.(S.func),'save') || ~isfield(S.(S.func).save,'suffix')
    S.(S.func).save.suffix = {''};
end

% get filename parts
S.(S.func).designmat = {'groups','subjects','sessions','blocks','conds'};
fnameparts = {'prefix','study','group','subject','session','block','cond','suffix','ext'};
parts = ismember(fnameparts,S.(S.func).fname.parts);

% get separators
partsep = {'_','_','_','_','_','_','_','_','_'};
if isfield(S.(S.func).fname,'partsep')
    partsep(parts) = S.(S.func).fname.partsep;
end

% load participant info and identify the relevant columns
[~,~,pdata] = xlsread(S.path.datfile);
grp_col = find(strcmp(pdata(1,:),'Group'));
sub_col = find(strcmp(pdata(1,:),'Subject'));
inc_col = find(strcmp(pdata(1,:),'Include'));

% identify number of groups of participants and find indices of
% participants in each group
if isnumeric([pdata{2:end,grp_col}])
    ugrp = unique([pdata{2:end,grp_col}]);
else
    grpdata = pdata(2:end,grp_col);
    grpdata = grpdata(~cell2mat(cellfun(@any,cellfun(@isnan,grpdata,'UniformOutput',false),'UniformOutput',false)));
    ugrp = unique(grpdata);
end
Ngrp = length(ugrp);
SubInd = cell(Ngrp,1);
Subs = [];
gn=0;
for g = 1:length(ugrp)
    gn=gn+1;
    inc_idx = find(cellfun(@(x) x>0, pdata(2:end,inc_col), 'UniformOutput', 1));
    if isnumeric([pdata{2:end,grp_col}])
        grp_idx = find(cellfun(@(x) x==ugrp(g), pdata(2:end,grp_col), 'UniformOutput', 1));
    else
        grp_idx = find(strcmp(pdata(2:end,grp_col),ugrp{g}));
    end
    SubInd{gn,1} = intersect(inc_idx,grp_idx);
    Nsub(gn,1) = length(SubInd{gn,1});
end

% make some corrections to the subject IDs if needed to make life easier later
if isnumeric(ugrp)
    ugrp = cellstr(num2str(ugrp'));
end
if isnumeric(S.(S.func).select.groups)
    S.(S.func).select.groups = cellfun(@num2str,S.(S.func).select.groups,'un',0);
end
grplist = find(ismember(ugrp,S.(S.func).select.groups)); %; % index of groups
if isempty(grplist)
    grplist = 1:Ngrp;
end
grps = ugrp(grplist);
SubInd=SubInd(grplist);
Ngrp=length(grps);

subjlists={};
for g = 1:Ngrp
    subgrp={};
    s2=0;
    for s = 1:Nsub(g)
        %if isnan(group(s,g))
        %    continue
        %end
        s2 = s2+1;
        subj = pdata{1+SubInd{g}(s),sub_col};
        if isnumeric(subj)
            subj = num2str(subj);
        end
        subgrp{s2,1} = subj;
    end
    subjlists{g,1} = subgrp;
end

%subjlists = subjlists(grplist);

% create file list
S.(S.func).filelist = {};
i = 0;
gs=0;
for g = 1:length(grps)
    for s = 1:length(subjlists{g,1}) 
        subj = subjlists{g,1}{s,1};
        subj_pdat_idx = SubInd{g}(s);
        grp = grps{g};
        gs=gs+1;
        
         % if subjects are specified in settings, only analyse those,
        % otherwise include all
        if ~isempty(S.(S.func).select.subjects{1}) 
            if ~any(strcmp(S.(S.func).select.subjects,subj))
                continue
            end
        else
            S.(S.func).subjects_in{gs} = subj;
        end
        
        if isempty(S.(S.func).study{:})
            study = {''};
        else
            study = {[S.(S.func).study{:} partsep{2}]};
        end
        if isempty(grp)
            grp = {''};
        else
            grp = {[grp partsep{3}]};
        end
        if isempty(subj)
            subj = {''};
        else
            subj = {[subj partsep{4}]};
        end
        if isempty(S.(S.func).load.prefix{:})
            pref = {''};
        else
            pref = {[S.(S.func).load.prefix{:} partsep{1}]};
        end
        if isempty(S.(S.func).load.suffix{:})
            suff = {''};
        else
            suff = S.(S.func).load.suffix;
        end
        if isempty(S.(S.func).fname.ext{:})
            ext = {'.*'};
        else
            ext = S.(S.func).fname.ext{:};
            if ~any(strfind(ext,'.'))
                ext = {['.' ext]};
            end
        end
        
        for a = 1:length(S.(S.func).select.sessions)
            for b = 1:length(S.(S.func).select.blocks)
                for c = 1:length(S.(S.func).select.conds)
                    
                    if isempty(S.(S.func).select.sessions{a})
                        session = {''};
                    else
                        session = {[S.(S.func).select.sessions{a} partsep{5}]};
                    end
                    if isempty(S.(S.func).select.blocks{b})
                        block = {''};
                    else
                        block = {[S.(S.func).select.blocks{b} partsep{6}]};
                    end
                    if isempty(S.(S.func).select.conds{c})
                        cond = {''};
                    else
                        cond = {[S.(S.func).select.conds{c} partsep{7}]};
                    end
                    
                    genname = [pref study grp subj session block cond suff ext];
                    genname = genname(parts);
                    genname = strjoin(genname,'');
                    
                    genname = strrep(genname,'*****','*');
                    genname = strrep(genname,'****','*');
                    genname = strrep(genname,'***','*');
                    genname = strrep(genname,'**','*');
                    genname = strrep(genname,'__.','.');
                    genname = strrep(genname,'_.','.');
                    %genname = regexprep(genname,{'\.','set'},{'','.set'});
                    
                    if isfield(S.(S.func),'grpdir')
                        grpdir = grps{g};
                    else
                        grpdir = '';
                    end
                    
                    if isfield(S.(S.func),'subjdir') && S.(S.func).subjdir
                        subdir = subjlists{g,1}{s,1};
                    else
                        subdir = '';
                    end
                    
                    disp(['loading: ' fullfile(S.path.file,grpdir,subdir,genname)])
                    file = dir(fullfile(S.path.file,grpdir,subdir,genname));
                    if length(file)~=1
                        disp(['No unique file named ' genname])
                        continue
                        %try 
                        %    file = file(1);
                        %    fname = file.name;
                        %catch
                        %    error('no files with this name')
                        %end
                    else
                        fname = file.name;
                    end
                    i=i+1;
                    S.(S.func).filelist{i} = fname;
                    S.(S.func).dirlist{i} = fullfile(S.path.file,grpdir,subdir);
                    S.(S.func).designmat{i+1,1} = grp{:}(1:end-1);
                    S.(S.func).designmat{i+1,2} = subj{:}(1:end-1);
                    S.(S.func).designmat{i+1,3} = S.(S.func).select.sessions{a};
                    S.(S.func).designmat{i+1,4} = S.(S.func).select.blocks{b};
                    S.(S.func).designmat{i+1,5} = S.(S.func).select.conds{c};
                    S.(S.func).subj_pdat_idx(i) = subj_pdat_idx;
                end
            end
        end
    end
end
S.(S.func).select.subjects = S.(S.func).subjects_in;

if isfield(S,'T')
    S=S.T;
end
