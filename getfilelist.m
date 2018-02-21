function S = getfilelist(S)

cd(S.filepath);

% if sessions/blocks/events are not specified in settings:
if isempty(S.sessions)
    S.sessions = {''};
end
if isempty(S.blocks)
    S.blocks = {''};
end
if isempty(S.conds)
    S.conds = {''};
end
if ~isfield(S,'loadprefix')
    S.loadprefix = '';
end
S.subjects_in = S.subjects;

% load participant info and identify the relevant columns
[~,~,pdata] = xlsread(S.datfile);
grp_col = find(strcmp(pdata(1,:),'Group'));
sub_col = find(strcmp(pdata(1,:),'Subject'));
inc_col = find(strcmp(pdata(1,:),'Include'));

% identify number of groups of participants and find indices of
% participants in each group
Ngrp = length(unique([pdata{2:end,grp_col}]));
SubInd = cell(Ngrp,1);
Subs = [];
gn=0;
for g = unique([pdata{2:end,grp_col}])
    gn=gn+1;
    inc_idx = find(cellfun(@(x) x==1, pdata(2:end,inc_col), 'UniformOutput', 1));
    grp_idx = find(cellfun(@(x) x==g, pdata(2:end,grp_col), 'UniformOutput', 1));
    SubInd{gn,1} = intersect(inc_idx,grp_idx);
    Nsub(gn,1) = length(SubInd{gn,1});
end

% make some corrections to the subject IDs if needed to make life easier later
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
grplist = 1:Ngrp; % index of groups

S.groups = subjlists(grplist);

% create file list
S.filelist = {};
i = 0;
gs=0;
for g = 1:length(S.groups)
    for s = 1:length(S.groups{g,1}) 
        subj = S.groups{g,1}{s,1};
        gs=gs+1;
        
        % if subjects are specified in settings, only analyse those,
        % otherwise include all
        if ~isempty(S.subjects) 
            if ~any(strcmp(S.subjects,subj))
                continue
            end
        else
            S.subjects_in{gs} = subj;
        end
        
        for a = 1:length(S.sessions)
            for b = 1:length(S.blocks)
                for c = 1:length(S.conds)
                    genname = [S.loadprefix '*' subj '*' S.sessions{a} '*' S.blocks{b} '*' S.conds{c} '*' S.loadext];
                    genname = strrep(genname,'**','*');
                    file = dir(fullfile(S.filepath,genname));
                    if length(file)~=1
                        error('File names are not uniquely specified')
                    else
                        fname = file.name;
                    end
                    i=i+1;
                    S.filelist{i} = fname;
                end
            end
        end
    end
end
S.subjects = S.subjects_in;
