%%--- function design_batch(D) ---%%

% D is an input structure requiring the following fields:
%--------------------------------------------------
%D.batch_path
%D.ffbatch
%D.data_path
%D.mask_path
%D.pdatfile
%D.spmstats_path
%D.anapref
%D.subdirpref
%D.subdirsuff
%D.imglist
%D.time_ana
%D.cond_list
%D.factors
%D.factortype
%D.interactions
%D.maineffects
%D.GMsca
%D.ancova
%D.grandmean
%D.globalnorm
%D.fcontrasts
%D.cov_names

% or for regression:
% D.regress_contrast
%%
%-----------------------------------------------------------------------  
function D=design_batch(D)
dbstop if error

if ~isfield(D,'anatype')
    D.anatype='group';
end
if ~isfield(D,'time_ana')
    D.time_ana=[];
end
if ~isfield(D,'identifier')
    D.identifier='';
end
if ~isfield(D,'baseline')
    D.baseline=[];
end
if ~isfield(D,'anapref')
    D.anapref='';
end
if ~isfield(D,'subdirsuff')
    D.subdirsuff='';
end
if ~isfield(D,'subdirpref')
    D.subdirpref='';
end
if ~isfield(D,'useIDfile')
    D.useIDfile=1;
end
if ~isfield(D,'askoverwrite')
    D.askoverwrite=1;
end
if isfield(D,'maineffects')
    isfactorial = 1;
else 
    isfactorial = 0;
    ismain = 0;
    isinter = 0;
end
if isfield(D,'regress_contrast')
    isregress = 1;
else
    isregress=0;
end
if ~isfield(D,'znorm')
    D.znorm = 0;
end
if ~isfield(D,'para')
    D.para = 0;
end
if ~isfield(D,'grandmean')
    D.grandmean = 0;
end
if ~isfield(D,'fileoptype')
    D.fileoptype = 'meancond';
end
if ~isfield(D,'timewin')
    D.timewin = 0;
end
if ~isfield(D,'meancentre')
    D.meancentre = 0;
end
if ~isfield(D,'pronto')
    D.pronto = 0;
end
if ~isfield(D,'imgpref')
    D.imgpref = '';
end
if ~isfield(D,'imgsuff')
    D.imgsuff = '';
end
if ~isfield(D,'maskfile')
    D.maskfile = '';
end

if D.para==1
    paraname = '_spm';
elseif D.para==2
    paraname = '_snpm';
elseif D.pronto
    paraname = '_prt';
end

% get factor names
factnames = '';
if isfactorial
    for f = 1:length(D.factors)
        factnames = [factnames '_' D.factors{f}];
    end
elseif isregress
    factnames = D.regress_contrastname;
end

covs = '';
for c = 1:length(D.cov_names)
    covs = [covs '_' D.cov_names{c}];
end
if ~isempty(D.time_ana)
    maskpref = '_m';
    for m = 1:length(D.time_ana)
        maskpref = [maskpref '_' num2str(D.time_ana(m))];
    end
else
    maskpref = '';
end
if iscell(D.identifier)
    identname = D.identifier{1};
else
    identname = D.identifier;
end
path_name = [D.anapref maskpref factnames covs D.subdirsuff paraname identname];
path_name = strrep(path_name,'*','');
D.spm_path = fullfile(D.spmstats_path,path_name); 
if exist(D.spm_path,'dir') && D.askoverwrite
    button = questdlg('Create new directory for analysis? ''No'' to overwrite');
    if strcmp(button,'Yes')
        dt = datestr(datetime,30);
        D.spm_path = [D.spm_path '_' dt];
    end
end
if ~exist(D.spm_path,'dir')
    mkdir(D.spm_path);
end
if isfield(D,'batch_path')
    copyfile(D.batch_path,D.spm_path);
end

if isfield(D,'TrainTest')
    TrainTestFac = D.factortype(find(~cellfun(@isempty,D.TrainTest)));
else
    TrainTestFac = {''};
end

[~,~,pdata] = xlsread(D.pdatfile);
if cell2mat(strfind(TrainTestFac,'w'))
    wlev = unique(D.TrainTest{ismember(D.factortype,'w')},'stable');
else
    wlev=0;
end
if cell2mat(strfind(TrainTestFac,'g'))
    glev = unique(D.TrainTest{ismember(D.factortype,'g')},'stable');
else
    glev=1;
end
if ~iscell(D.grphead); D.grphead = {D.grphead};end
if length(glev)>1 && length(D.grphead)>1
    grphead = D.grphead(glev);
else
    grphead = D.grphead(1);
end
    
for gh = 1:length(D.grphead)
    grp_col(gh) = find(strcmp(pdata(1,:),D.grphead{gh}));
end
sub_col = find(strcmp(pdata(1,:),D.subhead));
inc_col = find(strcmp(pdata(1,:),D.inchead));

if any(glev) && wlev
    % duplicate grouping to provide feature space for TrainTest
    grp_col = repmat(grp_col,1,2);
end

if isempty(grp_col) && ~isfield(D,'grpind')
    error(['no header in xls file called ' grphead]);
elseif isempty(sub_col)
    error(['no header in xls file called ' D.subhead]);
elseif isempty(inc_col)
    error(['no header in xls file called ' D.inchead]);
end

if isfield(D,'grpind')
    % grp indices already supplied
    grpdat = num2cell(D.grpind);
    grptype = unique([grpdat{:}]);
    Ngrp = length(grptype);
else
    for gd = 1:length(grp_col)
        grpdat = pdata(2:end,grp_col(gd));
        if isnumeric(grpdat{2,1})
            grpdat = grpdat(~cellfun(@isnan,grpdat));
            grptype = unique([grpdat{:}]);
            grptype = grptype(~isnan(grptype));
            grptype(grptype==0)=[];
            Ngrp(gd) = length(grptype);
        else %Change char Grp inputs to numbers
            grpdat = grpdat(~cellfun(@isnan,grpdat));
            grptype = unique(grpdat);
            grptype(isempty(grptype))=[];
            Ngrp(gd) = length(grptype);
            for g = 1:Ngrp(gd)
                grp_idx = cellfun(@(x) any(strcmp(grptype(g),x)), grpdat, 'UniformOutput', 0);
                grpdat(cell2mat(grp_idx)) = {[g]};
            end
        end
        grpdatall(:,gd) = grpdat;
    end
    grpdat = grpdatall;
end

% find index of subjects to include in analysis
inc_idx = cellfun(@(x) ismember(x,D.include_codes), pdata(2:end,inc_col), 'UniformOutput', 0);
inc_idx = find(cell2mat(inc_idx));

% find subject indices for each specific group
SubInd = {};
for gd = 1:length(Ngrp)
    for g = 1:Ngrp(gd)
        grp_idx = find(cellfun(@(x) x==g, grpdat(:,gd), 'UniformOutput', 1));
        if D.useIDfile==1
            SubInd{g,gd} = intersect(inc_idx,grp_idx);
        else
            SubInd{g,gd} = grp_idx;
        end
        Nsub(g,gd) = length(SubInd{g,gd});
    end
end

% select specific group(s)
if isfield(D,'grp_list')
   Nsub = Nsub(D.grp_list,:); 
   for gd = 1:length(Ngrp)
        Ngrp(gd) = length(Nsub(:,gd));
   end
   SubInd = SubInd(D.grp_list,:); 
end

ImgList={};
if isfield(D,'imglist_columnheaders')
    if ~isempty(D.imglist_columnheaders) && isempty(D.imglist)
        for ch = 1:length(D.imglist_columnheaders)
            ImgList(:,ch) = pdata(2:end,find(strcmp(pdata(1,:),D.imglist_columnheaders{ch})));
        end
        inum = cellfun(@isnumeric,ImgList);
        ImgList(inum) = cellfun(@num2str,ImgList(inum),'Uniformoutput',0);
    end
end

% select parametric or non-parametric analyses and the range of models
if isfactorial
    ismain = any(D.maineffects);
    isinter = any(D.interactions);
end
if D.para==1
    load(D.ffbatch);
elseif D.para==2
    load(D.npbatch);
    if isfactorial
        if ismain
            if strcmp(D.factortype(find(D.maineffects)),'w')
                matlabbatch = matlabbatch(1);
                model = 'w';
            elseif strcmp(D.factortype(find(D.maineffects)),'g')
                matlabbatch = matlabbatch(2);
                model = 'g';
            end
        elseif isinter
            uf = length(unique(D.factortype(find(D.interactions))));
            if uf==1 && strcmp(unique(D.factortype(find(D.interactions))),'w')% interaction of two within-subjects factors
                matlabbatch = matlabbatch(1); % paired test after subtraction of second factor
                model = 'ww';
            elseif uf==2 % interaction of within-subjects and group factor
                matlabbatch = matlabbatch(3); % unpaired test after subtraction of w factor
                model = 'gw';
            end
        else % assume one-sample t-test
            matlabbatch = matlabbatch(4);
            model = 'one';
        end
    end
elseif D.pronto
    load(D.batch);
    matlabbatch{1}.prt.data.group = [];
    paraname = '_pronto';
    if isfactorial
        if ismain
            if any(strcmp(D.factortype(find(D.maineffects)),'w')) && any(strcmp(D.factortype(find(D.maineffects)),'g')) && any(cell2mat(strfind(TrainTestFac,'g'))) && any(cell2mat(strfind(TrainTestFac,'w')))
                % training on within, testing on between, or vice versa
                %gn=0;
                %for gd = 1:length(Ngrp)
                    for g = 1:2 % 1=training "group", 2=testing "group"
                        %gn = gn+1;
                        matlabbatch{1}.prt.data.group(g).gr_name = ['Group' num2str(g)];
                    end
                %end
                model = 'wgtt';
            elseif strcmp(D.factortype(find(D.maineffects)),'w')
                if cell2mat(strfind(TrainTestFac,'g'))
                    for g = 1:Ngrp
                        matlabbatch{1}.prt.data.group(g).gr_name = ['Group' num2str(g)];
                    end
                    model = 'wg';
                else
                    matlabbatch{1}.prt.data.group(1).gr_name = 'all';
                    model = 'w';
                end
            elseif strcmp(D.factortype(find(D.maineffects)),'g')
                gn=0;
                model = '';
                for gd = 1:length(Ngrp)
                    for g = 1:Ngrp(gd)
                        gn = gn+1;
                        matlabbatch{1}.prt.data.group(gn).gr_name = ['Group' num2str(gn)];
                    end
                    model = [model 'g'];
                end
            end
        elseif isinter
            uf = length(unique(D.factortype(find(D.interactions))));
            if uf==1 && strcmp(unique(D.factortype(find(D.interactions))),'w')% interaction of two within-subjects factors
                matlabbatch{1}.prt.data.group(1).gr_name = 'all';
                model = 'ww';
            elseif uf==2 % interaction of within-subjects and group factor
                gn=0;
                model = 'w';
                for gd = 1:length(Ngrp)
                    for g = 1:Ngrp(gd)
                        gn = gn+1;
                        matlabbatch{1}.prt.data.group(gn).gr_name = ['Group' num2str(gn)];
                    end
                    model = ['g' model];
                end
            end
        end
    elseif isregress
        matlabbatch{1}.prt.data.group(1).gr_name = 'all';
    end
end


% add factors to the design

if D.para==1
    matlabbatch{1}.spm.stats.factorial_design.dir = {D.spm_path};
    if isfactorial
        for f = 1:length(D.factors)
            if strcmp(D.factortype(f),'w')
                dept = 1; %dependence
                variance = 0; % unequal variance
            elseif strcmp(D.factortype(f),'s')
                dept = 0; %dependence
                variance = 0; % unequal variance
            elseif strcmp(D.factortype(f),'g')
                dept = 0; %dependence
                variance = 1; % unequal variance
            end
            matlabbatch{1}.spm.stats.factorial_design.des.fblock.fac(f).name = D.factors{f};
            matlabbatch{1}.spm.stats.factorial_design.des.fblock.fac(f).dept = dept;
            matlabbatch{1}.spm.stats.factorial_design.des.fblock.fac(f).variance = variance;
            matlabbatch{1}.spm.stats.factorial_design.des.fblock.fac(f).gmsca = D.GMsca(f);
            matlabbatch{1}.spm.stats.factorial_design.des.fblock.fac(f).ancova = D.ancova(f);
        end
    end
elseif D.para==2
    generic.dir = {D.spm_path};
    generic.nPerm = D.nPerm;
    generic.vFWHM = D.vFWHM;
    generic.bVolm = D.bVolm;
    generic.ST.ST_U = D.ST_U;
    if isfactorial
        wfactind = double(ismember(D.factortype,'w'));
        if strcmp(model,'w')
            matlabbatch{1}.spm.tools.snpm.des.PairT = generic;
            scanname=D.factors(find(D.maineffects));
            scanname=scanname{:};
            wnum = find(wfactind .* D.maineffects);
        elseif strcmp(model,'ww')
            matlabbatch{1}.spm.tools.snpm.des.PairT = generic;
            scanname=horzcat(D.factors(find(D.interactions)));
            scanname=scanname{:};
            wnum = find(wfactind .* D.maineffects);
        elseif strcmp(model,'g')
            matlabbatch{1}.spm.tools.snpm.des.TwoSampT = generic;
            scanname='allavg';
            wnum = 0;
        elseif strcmp(model,'gw')
            matlabbatch{1}.spm.tools.snpm.des.TwoSampPairT = generic;
            wfactor = double(ismember(D.factortype,'w')) .* D.interactions;
            scanname=horzcat(D.factors(find(wfactor)));
            scanname=scanname{:};
            wnum = find(wfactind .* D.interactions);
        elseif strcmp(model,'one')
            matlabbatch{1}.spm.tools.snpm.des.OneSampT = generic;
            scanname='allavg';
            wnum = 0;
        end
    end
elseif D.pronto
    matlabbatch{1}.prt.data.dir_name = {D.spm_path};
    if isfactorial
        generic.mod_name = 'mod1';
        generic.TR = 0.0001;
        generic.design.new_design.unit = 0; % scans
        generic.design.new_design.multi_conds = {''};
        wfactind = double(ismember(D.factortype,'w'));
        gfactind = double(ismember(D.factortype,'g'));
        if strcmp(model,'w') || strcmp(model,'wg') || strcmp(model,'wgtt')
            wnum = find(wfactind .* D.maineffects);
            Nconds = length(unique(D.cond_list));
            for c = 1:Nconds
                generic.design.new_design.conds(c).cond_name = [D.factors{strcmp(D.factortype,'w')} num2str(c)];
                generic.design.new_design.conds(c).onsets = find(D.cond_list==c)'-1;
                generic.design.new_design.conds(c).durations = ones(length(generic.design.new_design.conds(c).onsets),1);
            end
        elseif strcmp(model,'g') || strcmp(model,'gg')
            wnum = 0;
            generic.design.new_design.conds(1).cond_name = 'allcond';
            generic.design.new_design.conds(1).onsets = find(D.cond_list)'-1;
            generic.design.new_design.conds(1).durations = ones(length(generic.design.new_design.conds(1).onsets),1);
        elseif strcmp(model,'ww')
            wnum = find(wfactind .* D.interactions);
            Nconds = length(unique(D.cond_list(:,wnum(1))));
            for c = 1:Nconds
                generic.design.new_design.conds(c).cond_name = [D.factors{wnum(1)} num2str(c)];
                generic.design.new_design.conds(c).onsets = c-1;
                generic.design.new_design.conds(c).durations = 1;%ones(length(generic.design.new_design.conds(c).onsets),1);
            end
            scanname=strcat(D.factors{wnum});
            %scanname=scanname{:};
        elseif strcmp(model,'gw') || strcmp(model,'ggw')
            wnum = find(wfactind .* D.interactions);
            gnum = find(gfactind .* D.interactions);
            generic.design.new_design.conds(1).cond_name = 'allcond';
            generic.design.new_design.conds(1).onsets = 0;%find(D.cond_list(:,1))'-1;
            generic.design.new_design.conds(1).durations = 1;%ones(length(generic.design.new_design.conds(1).onsets),1);
            scanname=strcat(D.factors{wnum});
            ucond = unique(D.cond_list);
            if length(ucond)==1
               scanname=[scanname '_'  num2str(ucond)];
            end
        end
    end 
end
if isregress
    scanname=D.regress_contrastname;
end
if D.znorm
    scanname_ext = '_znorm';
else
    scanname_ext = '';
end


% add main effects and interactions
if isfactorial
    if D.para==1
        ne = 0;
        for e = 1:size(D.interactions,1)
            ne = ne+1;
            matlabbatch{1}.spm.stats.factorial_design.des.fblock.maininters{ne}.inter.fnums = find(D.interactions)';
        end
        me_ind = find(D.maineffects);
        for e = 1:length(me_ind)
            ne = ne+1;
            matlabbatch{1}.spm.stats.factorial_design.des.fblock.maininters{ne}.fmain.fnum = me_ind(e);
        end
    end
end


% find and load scans into matlabbatch
gs = 0;
subID={};
for gd = 1:length(Ngrp)
    gds = 0;
    for g = 1:Ngrp(gd)
        for s = 1:Nsub(g,gd)
            gs = gs+1;
            gds = gds+1;
            if D.useIDfile==1
                subID{gs} = deblank(num2str((pdata{SubInd{g,gd}(s)+1,sub_col})));
                if isnumeric(subID{gs}); subID{gs} = num2str(subID{gs}); end;
                if D.folder
                    subdir = fullfile(D.data_path, [D.anapref D.subdirpref subID{gs} D.subdirsuff]);
                    if ~exist(subdir,'dir')
                        subdir = dir(fullfile(D.data_path, [D.anapref D.subdirpref subID{gs} D.subdirsuff]));
                        if length(subdir)>1
                            error('Subject directory name is not unique in this folder');
                            display(fullfile(D.data_path, [D.anapref D.subdirpref subID{gs} D.subdirsuff]));
                        elseif length(subdir)==0
                            error('no directory found with this name');
                            display(fullfile(D.data_path, [D.anapref D.subdirpref subID{gs} D.subdirsuff]));
                        end
                        subdir=fullfile(D.data_path, subdir.name);
                    end
                    subfile = '';
                else
                    subdir = D.data_path;
                    subfile = [D.anapref D.subdirpref subID{gs} D.subdirsuff];
                end
            else
                subdir = D.data_path;
                subfile = D.imglist(gs,:)';
            end
            
            if ~isempty(ImgList)
                D.imglist = ImgList(SubInd{g,gd}(s),:)';
            end
            if ~isempty(D.imgpref)
                D.imglist = strcat(D.imgpref,D.imglist);
            end
            if ~isempty(D.imgsuff)
                D.imglist = strcat(D.imglist,D.imgsuff);
            end

            subimg = cell(1,1);
            subimgcond=[];
            if isfactorial
                if D.para==1 || (D.pronto && ~isinter)
                    % for spm, load all files
                    ii=0;
                    for i = 1:length(D.imglist)
                        fname = {fullfile(subdir, [subfile D.imglist{i}])};
                        if ~exist(fname{:},'file')
                            fname_temp = dir(fullfile(subdir, [subfile D.imglist{i}]));
                            if length(fname_temp)>1 && ~strcmp(D.anatype,'singlesubject')
                                display(fname);
                                error('Subject file name is not unique in this folder');
                            elseif length(fname_temp)==0
                                display(fname);
                                error('no file found with this name');
                            end
                            for f = 1:length(fname_temp)
                                fname{f}=fullfile(subdir, fname_temp(f).name);
                            end
                        end
                        for f = 1:length(fname)
                            ii = ii+1;
                            subimg{ii,1} = [fname{f} ',1'];
                            subimgcond(ii,1) = i;
                        end
                    end
                    if strcmp(D.anatype,'singlesubject')
                        condind=[];
                        for c = 1:Nconds
                            ntrial(c) = [condind sum(ismember(subimgcond,find(D.cond_list==c)))];
                            generic.design.new_design.conds(c).durations = ones(ntrial(c),1);
                            generic.design.new_design.conds(c).onsets = find(ismember(subimgcond,find(D.cond_list==c)))'-1;
                            generic.design.new_design.conds(c).blocks = generic.design.new_design.conds(c).onsets+1;
                        end
                        %if D.randomise
                        %    blockind = randperm(sum(ntrial));
                        %    condind = [0 cumsum(ntrial)];
                        %    subimg = subimg(blockind);
                        %    for c = 1:Nconds
                        %        cind = condind(c)+1:condind(c+1);
                        %        generic.design.new_design.conds(c).onsets = blockind(ismember(blockind,cind))-1;
                        %        %generic.design.new_design.conds(c).blocks = blockind(ismember(blockind,cind));
                        %    end
                        %end
                    end
                elseif D.para==2 || (D.pronto && isinter)
                    if D.useIDfile==1
                        % for SnPM, check averaged files exist; if not then create them.
                        fnames = dir(fullfile(subdir, [scanname '*' D.fileoptype '*' scanname_ext '*']));
                        if wnum==0
                            condlist = ones(length(D.cond_list),1);
                        else
                            condlist = D.cond_list(:,find(find(wfactind)==wnum));
                        end
                        if isempty(fnames) || D.overwrite
                            factor_img(subdir,D.imglist,condlist,scanname,D.fileoptype,D.znorm);
                            fnames = dir(fullfile(subdir, [scanname '*' D.fileoptype '*' scanname_ext '*']));
                        end
                        for i=1:length(fnames)
                            subimg{i,1} = fullfile(subdir,[fnames(i).name ',1']);
                        end
                    else
                        subimg = subfile;
                    end
                end

                % if there is a Group factor, add in group level to cond_list
                grpfactind = find(strcmp(D.factortype,'g'));
                if ~isempty(grpfactind)
                    nlist = D.cond_list;
                    if grpfactind<=size(D.cond_list,2)
                        nlist(:,grpfactind+1:end+1) = nlist(:,grpfactind:end);
                    end
                    nlist(:,grpfactind) = g*ones(size(nlist,1),1);
                    condlist = nlist;
                else
                    condlist = D.cond_list;
                end

            elseif isregress
                if D.useIDfile==1
                    % check averaged files exist; if not then create them.
                    %fnames = dir(fullfile(subdir, ['*' scanname scanname_ext '*']));
                    condlist = D.regress_contrast;
                    if length(D.imglist)>1
                        fname = factor_img(subdir,D.imglist,condlist,scanname,D.fileoptype,D.znorm);
                        %fnames = dir(fullfile(subdir, ['*' scanname scanname_ext '*']));
                    else
                        fname = dir(fullfile(subdir,D.imglist{:}));
                        fname = fullfile(subdir,fname.name);
                    end
                    for i=1%:length(fnames)
                        subimg{i,1} = [fname ',1'];
                    end
                else
                    subimg = subfile;
                end
            end
            scans = subimg;

            % add the above scan and group info to the matlabbatch, or create new scans
            if isfactorial
                if D.para==1
                    matlabbatch{1}.spm.stats.factorial_design.des.fblock.fsuball.fsubject(gs).scans = scans;
                    matlabbatch{1}.spm.stats.factorial_design.des.fblock.fsuball.fsubject(gs).conds = condlist;
                elseif D.para==2
                    if strcmp(model,'w')
                        matlabbatch{1}.spm.tools.snpm.des.PairT.fsubject(gs).scans = scans;
                        matlabbatch{1}.spm.tools.snpm.des.PairT.fsubject(gs).scindex = unique(condlist,'stable');
                    elseif strcmp(model,'ww')
                        %matlabbatch{1}.spm.tools.snpm.des.PairT.fsubject(gs).scans = ;
                        %matlabbatch{1}.spm.tools.snpm.des.PairT.fsubject(gs).scindex = ;
                    elseif strcmp(model,'g')
                        if g==1
                            matlabbatch{1}.spm.tools.snpm.des.TwoSampT.scans1(s,1) = scans;
                        elseif g==2
                            matlabbatch{1}.spm.tools.snpm.des.TwoSampT.scans2(s,1) = scans;
                        end
                    elseif strcmp(model,'gw')
                        if g==1
                            %matlabbatch{1}.spm.tools.snpm.des.TwoSampPairT.scans1.fsubject(s).scans = ;
                            %matlabbatch{1}.spm.tools.snpm.des.TwoSampPairT.scans1.fsubject(s).scindex = ;
                        elseif g==2
                            %matlabbatch{1}.spm.tools.snpm.des.TwoSampPairT.scans2.fsubject(s).scans = ;
                            %matlabbatch{1}.spm.tools.snpm.des.TwoSampPairT.scans2.fsubject(s).scindex = ;
                        end
                    elseif strcmp(model,'one')
                        matlabbatch{1}.spm.tools.snpm.des.OneSampT.P(gs,1) = scans;
                    end
                elseif D.pronto
                    if strcmp(model,'wgtt')
                        if gd==1 % within-subject factor "group"
                            matlabbatch{1}.prt.data.group(1).select.subject{gds,1} = generic;
                            matlabbatch{1}.prt.data.group(1).select.subject{gds,1}.scans = scans;
                        elseif gd==2 % between-subject factor "group"
                            %generic.design.new_design.conds(:) = generic.design.new_design.conds(D.grpcond);
                            %generic.design.new_design.conds.cond_name = [D.factors{wnum(1)} num2str(g)]; % change the name to match the within levels so that pronto can perform testing
                            scans = [scans(D.cond_list==D.grpcond);scans(D.cond_list==D.grpcond)];
                            matlabbatch{1}.prt.data.group(2).select.subject{gds,1} = generic;
                            matlabbatch{1}.prt.data.group(2).select.subject{gds,1}.scans = scans;
                        end
                    elseif strcmp(model,'g') || strcmp(model,'gw') || strcmp(model,'gg') || strcmp(model,'ggw') || strcmp(model,'wg')
                         matlabbatch{1}.prt.data.group(g+2*(gd-1)).select.subject{s,1} = generic;
                         matlabbatch{1}.prt.data.group(g+2*(gd-1)).select.subject{s,1}.scans = scans;
                    elseif strcmp(model,'w') || strcmp(model,'ww')
                        matlabbatch{1}.prt.data.group(1).select.subject{gs,1} = generic;
                        matlabbatch{1}.prt.data.group(1).select.subject{gs,1}.scans = scans;
                    end
                end
            elseif isregress
                if D.para==1
                    matlabbatch{1}.spm.stats.factorial_design.des.mreg.scans(gs,1) = scans;
                elseif D.para==2
                    matlabbatch{1}.spm.tools.snpm.des.Corr.P(gs,1) = scans;
                elseif D.pronto
                    matlabbatch{1}.prt.data.group(g).select.modality.mod_name = 'mod1';
                    matlabbatch{1}.prt.data.group(g).select.modality.subjects(s,1) = scans;
                end
            end
        end
    end
end

% covariates
if ~isempty(D.cov_names)
    if ~strcmp(D.cov_names{1},'')
        for c = 1:length(D.cov_names)
            cov_col = find(strcmp(pdata(1,:),D.cov_names{c}));
            covdat = pdata(2:end,cov_col);

            % convert categorical (char) covariate data to numbers
            if all(cellfun(@ischar,covdat))
                covtype = unique(covdat);
                Ncovtype = length(covtype);
                for c = 1:Ncovtype
                    cov_idx = cellfun(@(x) any(strcmp(covtype(c),x)), covdat, 'UniformOutput', 0);
                    covdat(cell2mat(cov_idx)) = {[c]};
                end
            elseif any(cellfun(@ischar,covdat))
                covdat(cellfun(@ischar,covdat)) = repmat({NaN},sum(cellfun(@ischar,covdat)),1);
            end

            covdat = [covdat{:}]';
            covdat = covdat(vertcat(SubInd{:}));
            covdat = repmat(covdat,1,length(scans));
            covdat = reshape(covdat', size(covdat,1)*size(covdat,2),1);
            covdat(isnan(covdat)) = mean(covdat(~isnan(covdat)));
            if isfactorial
                if D.para==1
                    matlabbatch{1}.spm.stats.factorial_design.cov(c).c = covdat; % vector
                    matlabbatch{1}.spm.stats.factorial_design.cov(c).cname = D.cov_names{c}; % name
                    matlabbatch{1}.spm.stats.factorial_design.cov(c).iCFI = 1;% interaction with factor
                    matlabbatch{1}.spm.stats.factorial_design.cov(c).iCC = 1;% centering
                elseif D.para==2
                    if strcmp(model,'g')
                        matlabbatch{1}.spm.tools.snpm.des.TwoSampT.cov(c).c = covdat; % vector
                        matlabbatch{1}.spm.tools.snpm.des.TwoSampT.cov(c).cname = D.cov_names{c}; % name
                    elseif strcmp(model,'gw')
                        matlabbatch{1}.spm.tools.snpm.des.TwoSampPairT.cov(c).c = covdat; % vector
                        matlabbatch{1}.spm.tools.snpm.des.TwoSampPairT.cov(c).cname = D.cov_names{c}; % name
                    elseif strcmp(model,'one')
                        matlabbatch{1}.spm.tools.snpm.des.OneSampT.cov(c).c = covdat; % vector
                        matlabbatch{1}.spm.tools.snpm.des.OneSampT.cov(c).cname = D.cov_names{c}; % name
                    end
                end
            elseif isregress
                if D.para==1
                    matlabbatch{1}.spm.stats.factorial_design.des.mreg.mcov(c).c = covdat; % vector
                    matlabbatch{1}.spm.stats.factorial_design.des.mreg.mcov(c).cname = D.cov_names{c}; % name
                elseif D.para==2
                    matlabbatch{1}.spm.tools.snpm.des.Corr.CovInt = covdat; % vector
                elseif D.pronto
                    gs = 0;
                    for g = 1:Ngrp
                        for s = 1:Nsub(g)
                            gs = gs+1;
                            matlabbatch{1}.prt.data.group(g).select.modality.rt_subj(s,1) = covdat(gs); % vector
                        end
                    end
                end
            end
        end
    end
else
    if D.para==1
        matlabbatch{1}.spm.stats.factorial_design.cov = struct('c', {}, 'cname', {}, 'iCFI', {}, 'iCC', {}); 
    elseif D.para==2
        if isfactorial
            if strcmp(model,'g')
                matlabbatch{1}.spm.tools.snpm.des.TwoSampT.cov = struct('c', {}, 'cname', {}); 
            elseif strcmp(model,'gw')
                matlabbatch{1}.spm.tools.snpm.des.TwoSampPairT.cov = struct('c', {}, 'cname', {});
            end
        elseif isregress
            matlabbatch{1}.spm.tools.snpm.des.Corr.cov = struct('c', {}, 'cname', {});
        end
    end
end
if isfactorial
    if D.para==1
        matlabbatch{1}.spm.stats.factorial_design.multi_cov = struct('files', {}, 'iCFI', {}, 'iCC', {});
    end
end

% specify mask
if ~isempty(D.time_ana)
    maskname = [num2str(D.time_ana(1)) '_' num2str(D.time_ana(2)) '.img'];
    maskfile = fullfile(D.mask_path,maskname);
    %if ~exist(maskfile,'file')
        S = struct;
        S.image = subimg{1}(1:end-2); % example image
        S.timewin = D.time_ana;
        S.outfile = maskfile;
        spm_eeg_mask(S)
    %end
    masking.em = {[maskfile  ',1']};
elseif ~isempty(D.maskfile)
    masking.em = {[D.maskfile  ',1']};
else
    masking.em = {''};
end

% specify ROI atlas
if D.pronto
    ROIfile='';
    if ~isempty(D.timewin) && ~isempty(D.time_ana)
        % create window range
        winr=[];
        for i = 1:Inf
            if D.time_ana(1)+D.timewin(1)-1+(D.timewin(1)*(i-1))>D.time_ana(2)
                break
            else
                winr(i,:) = [D.time_ana(1),D.time_ana(1)+D.timewin(1)-1]+(D.timewin(1)*(i-1));
            end
        end
        ROIname = [num2str(D.time_ana(1)) '_' num2str(D.time_ana(2)) '_step' num2str(D.timewin) '.img'];
        ROIfile = fullfile(D.mask_path,ROIname);
        V=spm_vol(maskfile);
        Y = spm_read_vols(V);
        Nt=size(Y,3);
        for i = 1:size(winr,1)
            begsample = V.mat\[0 0 winr(i,1) 1]';
            begsample = begsample(3);
            endsample = V.mat\[0 0 winr(i,2) 1]';
            endsample = endsample(3);
            if any([begsample endsample] < 0) || any([begsample endsample] > Nt)
                error('The window is out of limits for the image.');
            end
            [junk,begsample] = min(abs(begsample-(1:Nt)));
            [junk,endsample] = min(abs(endsample-(1:Nt)));
            Y(: , :, begsample:endsample)  = i;
            if i==1
                Y(: , :, 1:(begsample-1))   = 0;
            elseif i==size(winr,1)
                Y(: , :, (endsample+1):end) = 0;
            end
        end
        V.fname=ROIfile;
        spm_write_vol(V,Y);
    elseif isfield(D,'ROIfile')
        ROIfile = D.ROIfile;
    end
end

masking.tm.tm_none = 1;
masking.im = 1;
if D.para==1
    matlabbatch{1}.spm.stats.factorial_design.masking = masking;
    matlabbatch{1}.spm.stats.factorial_design.globalc.g_omit = 1;
    %globalm=matlabbatch{1}.spm.stats.factorial_design.globalm;
elseif D.para==2
    if isfactorial
        if strcmp(model,'w') || strcmp(model,'ww')
            matlabbatch{1}.spm.tools.snpm.des.PairT.masking = masking;
            matlabbatch{1}.spm.tools.snpm.des.PairT.globalc.g_omit = 1;
            %globalm=matlabbatch{1}.spm.tools.snpm.des.PairT.globalm;
        elseif strcmp(model,'g')
            matlabbatch{1}.spm.tools.snpm.des.TwoSampT.masking = masking;
            matlabbatch{1}.spm.tools.snpm.des.TwoSampT.globalc.g_omit = 1;
            %globalm=matlabbatch{1}.spm.tools.snpm.des.TwoSampT.globalm;
        elseif strcmp(model,'gw')
            matlabbatch{1}.spm.tools.snpm.des.TwoSampPairT.masking = masking;
            matlabbatch{1}.spm.tools.snpm.des.TwoSampPairT.globalc.g_omit = 1;
            %globalm=matlabbatch{1}.spm.tools.snpm.des.TwoSampPairT.globalm;
        elseif strcmp(model,'one')
            matlabbatch{1}.spm.tools.snpm.des.OneSampPairT.masking = masking;
            matlabbatch{1}.spm.tools.snpm.des.OneSampPairT.globalc.g_omit = 1;
            %globalm=matlabbatch{1}.spm.tools.snpm.des.TwoSampPairT.globalm;
        end
    elseif isregress
        matlabbatch{1}.spm.tools.snpm.des.Corr.masking = masking;
        matlabbatch{1}.spm.tools.snpm.des.Corr.globalc.g_omit = 1;
    end
elseif D.pronto
    matlabbatch{1}.prt.data.mask.mod_name = 'mod1';
    matlabbatch{1}.prt.data.mask.fmask = masking.em;
end


%specify scaling parameters
if D.para
    if D.grandmean
        globalm.gmsca.gmsca_yes.gmscv = D.grandmean;
        %if isfield(globalm.gmsca, 'gmsca_no')
        %    globalm.gmsca = rmfield(globalm.gmsca, 'gmsca_no');
        %end
    else
        globalm.gmsca.gmsca_no = 1;
        %if isfield(globalm.gmsca, 'gmsca_yes')
        %    globalm.gmsca = rmfield(globalm.gmsca, 'gmsca_yes');
        %end
    end
    globalm.glonorm = D.globalnorm;
end

if D.para==1
    matlabbatch{1}.spm.stats.factorial_design.globalm = globalm;

    matlabbatch{2}.spm.stats.fmri_est.spmmat(1) = cfg_dep('Factorial design specification: SPM.mat File', substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));
    matlabbatch{2}.spm.stats.fmri_est.write_residuals = D.resid;
    matlabbatch{2}.spm.stats.fmri_est.method.Classical = 1;

    matlabbatch{3} = struct();
    matlabbatch{3}.spm.stats.con.spmmat = {fullfile(D.spm_path,'SPM.mat')};
    matlabbatch{3}.spm.stats.con.delete = 1;
    matlabbatch{3}.spm.stats.con.consess = cell(1,size(D.fcontrasts,1)+size(D.tcontrasts,1));
    for fc = 1:size(D.fcontrasts,1)
        matlabbatch{3}.spm.stats.con.consess{1,fc}.fcon.name = D.fcontrasts{fc,2};
        matlabbatch{3}.spm.stats.con.consess{1,fc}.fcon.weights = D.fcontrasts{fc,1};
        matlabbatch{3}.spm.stats.con.consess{1,fc}.fcon.sessrep = 'none';
    end
    tcn=0;
    for tc = fc+1:size(D.fcontrasts,1)+size(D.tcontrasts,1)
        tcn=tcn+1;
        matlabbatch{3}.spm.stats.con.consess{1,tc}.tcon.name = D.tcontrasts{tcn,2};
        matlabbatch{3}.spm.stats.con.consess{1,tc}.tcon.weights = D.tcontrasts{tcn,1};
        matlabbatch{3}.spm.stats.con.consess{1,tc}.tcon.sessrep = 'none';
    end
elseif D.para==2
    if isfactorial
        if strcmp(model,'w') || strcmp(model,'ww')
            matlabbatch{1}.spm.tools.snpm.des.PairT.globalm = globalm;
        elseif strcmp(model,'g')
            matlabbatch{1}.spm.tools.snpm.des.TwoSampT.globalm = globalm;
        elseif strcmp(model,'gw')
            matlabbatch{1}.spm.tools.snpm.des.TwoSampPairT.globalm = globalm;
        elseif strcmp(model,'one')
            matlabbatch{1}.spm.tools.snpm.des.OneSampT.globalm = globalm;
        end
        if isregress
            matlabbatch{1}.spm.tools.snpm.des.Corr.globalm = globalm;
        end
    end

    matlabbatch{2}.spm.tools.snpm.cp.snpmcfg = {fullfile(D.spm_path,'SnPMcfg.mat')};

    matlabbatch{3} = struct();
    matlabbatch{3}.spm.tools.snpm.inference.SnPMmat = {fullfile(D.spm_path,'SnPM.mat')};
    matlabbatch{3}.spm.tools.snpm.inference.Thr.Clus.ClusSize.CFth = NaN; % set to NaN if already defined in the design
    matlabbatch{3}.spm.tools.snpm.inference.Thr.Clus.ClusSize.ClusSig.FWEthC = 0.05;
    matlabbatch{3}.spm.tools.snpm.inference.Tsign=1;
    matlabbatch{3}.spm.tools.snpm.inference.WriteFiltImg.name='SnPM_filtered';
    matlabbatch{3}.spm.tools.snpm.inference.Report = 'MIPtable';
elseif D.pronto
    % kernel / feature selection
    matlabbatch{2}.prt.fs.k_file = 'kernelname';
    matlabbatch{2}.prt.fs.modality.voxels=[];
    matlabbatch{2}.prt.fs.modality.voxels.all_voxels = 1;
    matlabbatch{2}.prt.fs.modality.atlasroi{1, 1} = '';
    if ~isempty(ROIfile)
        matlabbatch{2}.prt.fs.modality.atlasroi{1,1}   = [ROIfile ',1'];
    end
    matlabbatch{2}.prt.fs.infile = {fullfile(D.spm_path,'PRT.mat')};
    matlabbatch{2}.prt.fs.modality.mod_name = 'mod1';
    matlabbatch{2}.prt.fs.modality.conditions = [];
    if isfactorial
        matlabbatch{2}.prt.fs.modality.conditions.all_cond = 1;
    elseif isregress
        matlabbatch{2}.prt.fs.modality.conditions.all_scans = 1;
    end
    matlabbatch{2}.prt.fs.modality.detrend.no_dt = 1; % no detrend
    matlabbatch{2}.prt.fs.modality.normalise.no_gms = 1;
    matlabbatch{2}.prt.fs.flag_mm = 0;
    
    
    % model
    matlabbatch{3}.prt.model.infile = {fullfile(D.spm_path,'PRT.mat')};
    matlabbatch{3}.prt.model.use_kernel = D.kernel;
    matlabbatch{3}.prt.model.fsets = 'kernelname';
    matlabbatch{3}.prt.model.model_type = [];
    matlabbatch{3}.prt.model.cv_type = [];
    CV=[];
    if isfactorial
        try
            matlabbatch{3}.prt.model.model_name = D.factors{find(D.maineffects)};
        catch
            matlabbatch{3}.prt.model.model_name = D.factors{find(D.interactions(1))};
        end
        matlabbatch{3}.prt.model.model_type.classification.machine_cl = [];
        switch D.machine 
            case 'svm_binary'
                matlabbatch{3}.prt.model.model_type.classification.machine_cl.svm.svm_opt = 1; % optimise hyperparameter?
                matlabbatch{3}.prt.model.model_type.classification.machine_cl.svm.svm_args = [0.01 0.1 1 10 100]; % soft-margin hyperparameter range (only if svm_opt is 1)
                matlabbatch{3}.prt.model.model_type.classification.machine_cl.svm.cv_type_nested = [];
                matlabbatch{3}.prt.model.model_type.classification.machine_cl.svm.cv_type_nested.(D.cv_type) = 1;
            case 'gpc_binary'
                matlabbatch{3}.prt.model.model_type.classification.machine_cl.gpc.gpc_args = '-l erf -h';
            case 'gpc_multi'
                error('PRT_MODEL.M HAS BEEN MODIFIED TO ONLY ALLOW BINARY')
                matlabbatch{3}.prt.model.model_type.classification.machine_cl.gpclap.gpclap_args = '-h';
            case 'mkl'
                matlabbatch{3}.prt.model.model_type.classification.machine_cl.sMKL_cla.sMKL_cla_opt = 1;
                matlabbatch{3}.prt.model.model_type.classification.machine_cl.sMKL_cla.sMKL_cla_args = [0.01 0.1 1 10 100];
                matlabbatch{3}.prt.model.model_type.classification.machine_cl.sMKL_cla.cv_type_nested = [];
                matlabbatch{3}.prt.model.model_type.classification.machine_cl.sMKL_cla.cv_type_nested.(D.cv_type) = 1;
        end
        
        matlabbatch{3}.prt.model.include_allscans = 0;
        if strcmp(model,'w') || strcmp(model,'ww')
            % classes
            Nconds = length(unique(D.cond_list(:,1)));
            for c = 1:Nconds
                matlabbatch{3}.prt.model.model_type.classification.class(c).class_name = [D.factors{wnum} num2str(c)];
                matlabbatch{3}.prt.model.model_type.classification.class(c).group.gr_name = 'all';
                matlabbatch{3}.prt.model.model_type.classification.class(c).group.subj_nums = (1:sum(Nsub))';
                matlabbatch{3}.prt.model.model_type.classification.class(c).group.conditions=[];
                matlabbatch{3}.prt.model.model_type.classification.class(c).group.conditions.conds.cond_name = [D.factors{wnum(1)} num2str(c)];
            end
        elseif strcmp(model,'g') || strcmp(model,'gw') || strcmp(model,'gg') || strcmp(model,'ggw')
            % classes
            c=0;
            for gd = 1:length(Ngrp)
                for g = 1:Ngrp(gd)
                    c=c+1;
                    matlabbatch{3}.prt.model.model_type.classification.class(c).class_name = ['Group' num2str(c)];
                    matlabbatch{3}.prt.model.model_type.classification.class(c).group.gr_name = ['Group' num2str(c)];
                    matlabbatch{3}.prt.model.model_type.classification.class(c).group.subj_nums = (1:sum(Nsub(g,gd)))'; %(subind(c)+1:subind(c+1))';
                    matlabbatch{3}.prt.model.model_type.classification.class(c).group.conditions=[];
                    matlabbatch{3}.prt.model.model_type.classification.class(c).group.conditions.all_cond = 1;
                end
            end
        elseif strcmp(model,'wg')
            % classes
            Nconds = length(unique(D.cond_list));
            for g = 1:Ngrp
                for c = 1:Nconds
                    matlabbatch{3}.prt.model.model_type.classification.class(c).class_name = [D.factors{wnum} num2str(c)];
                    matlabbatch{3}.prt.model.model_type.classification.class(c).group(g).gr_name = ['Group' num2str(g)];
                    matlabbatch{3}.prt.model.model_type.classification.class(c).group(g).subj_nums = (1:Nsub(g))';
                    matlabbatch{3}.prt.model.model_type.classification.class(c).group(g).conditions=[];
                    matlabbatch{3}.prt.model.model_type.classification.class(c).group(g).conditions.conds.cond_name = [D.factors{wnum} num2str(c)];
                end
            end
        elseif strcmp(model,'wgtt')
            Nconds = length(unique(D.cond_list));
            %if wlev==1 && glev==2
                % class indices
           %     cw = 1:Nconds;
           %     gw = Nconds+[1:Ngrp(2)];
           % elseif wlev==2 && glev==1
           %     gw = 1:Ngrp(1);
           %     cw = Ngrp(1)+[1:Nconds];
           % end
            % within factor
            for c = 1:Nconds
                matlabbatch{3}.prt.model.model_type.classification.class(c).class_name = [D.factors{wnum} num2str(c)];
                matlabbatch{3}.prt.model.model_type.classification.class(c).group.gr_name = 'Group1';
                matlabbatch{3}.prt.model.model_type.classification.class(c).group.subj_nums = (1:sum(Nsub(:,1)))';
                matlabbatch{3}.prt.model.model_type.classification.class(c).group.conditions=[];
                matlabbatch{3}.prt.model.model_type.classification.class(c).group.conditions.conds.cond_name = [D.factors{wnum} num2str(c)];
                
            end
            % between factor
            subind = [0;cumsum(Nsub(:,1))];
            for c = 1:length(Ngrp)
                cn=c+Nconds;
                matlabbatch{3}.prt.model.model_type.classification.class(cn).class_name = ['Group' num2str(c)];
                matlabbatch{3}.prt.model.model_type.classification.class(cn).group.gr_name = 'Group2';
                matlabbatch{3}.prt.model.model_type.classification.class(cn).group.subj_nums = (subind(c)+1:subind(c+1))';
                matlabbatch{3}.prt.model.model_type.classification.class(cn).group.conditions=[];
                if D.grpcond
                %    matlabbatch{3}.prt.model.model_type.classification.class(cn).group.conditions.conds.cond_name = [D.factors{wnum} num2str(c)];
                %else
                    matlabbatch{3}.prt.model.model_type.classification.class(cn).group.conditions.all_cond = 1;
                %    CV = [CV; glev*ones(Nsub(g,2)*Nconds,1)];
                end
            end
            %CV = [CV; wlev*ones(sum(Nsub(:,1))*length(find(D.cond_list)),1)];
            %CV = [CV; glev*ones(Nsub(g,2)*length(find(D.cond_list==D.grpcond)),1)];
        end
        
        % factorial CV
        if cell2mat(strfind(TrainTestFac,'g'))
            if isempty(CV)
                if length(Ngrp)>1 % train and test on different groupings, or across within and between factors
                    for gd = 1:length(Ngrp)
                        CV = [CV;gd*ones(sum(Nsub(:,gd))*length(scans),1)];
                    end
                elseif length(Ngrp)==1 % train and test on different levels within a single grouping
                    for g = 1:Ngrp
                        CV = [CV;glev(g)*ones(Nsub(g)*length(scans),1)];
                    end
                end
            end
            save(fullfile(D.spm_path,'custom_CV.mat'),'CV');
        end
    
    elseif isregress
        matlabbatch{3}.prt.model.model_name = D.regress_contrastname;
        matlabbatch{3}.prt.model.model_type.regression.machine_rg = [];
        switch D.machine 
            case 'gpr'
                matlabbatch{3}.prt.model.model_type.regression.machine_rg.gpr.gpr_args = '-l gauss -h';
            case 'krr'
                matlabbatch{3}.prt.model.model_type.regression.machine_rg.krr.krr_opt = 1;
                matlabbatch{3}.prt.model.model_type.regression.machine_rg.krr.krr_args = 1;
                matlabbatch{3}.prt.model.model_type.regression.machine_rg.krr.cv_type_nested = [];
                matlabbatch{3}.prt.model.model_type.regression.machine_rg.krr.cv_type_nested.(D.cv_type) = 1;
        end
       
        matlabbatch{3}.prt.model.include_allscans = 0;
        
        % groups
        subind = [0;cumsum(Nsub)];
        for c = 1:Ngrp
            matlabbatch{1}.prt.data.group(c).gr_name = ['Group' num2str(c)];
            matlabbatch{3}.prt.model.model_type.regression.reg_group(c).gr_name = ['Group' num2str(c)];
            matlabbatch{3}.prt.model.model_type.regression.reg_group(c).subj_nums = (1:sum(Nsub(c)))'; %(subind(c)+1:subind(c+1))';
        end
    end
    
    % cross-validation
    if ~isempty(CV)
        matlabbatch{3}.prt.model.cv_type.cv_custom = {fullfile(D.spm_path,'custom_CV.mat')}; 
    elseif strcmp(D.cv_type,'cv_lkso')
        matlabbatch{3}.prt.model.cv_type.cv_lkso.k_args = D.nfolds;
    elseif strcmp(D.cv_type,'cv_lkbo')
        matlabbatch{3}.prt.model.cv_type.cv_lkbo.k_args = D.nfolds;
    else
        matlabbatch{3}.prt.model.cv_type.(D.cv_type) = 1;
    end
    
    % operations
    matlabbatch{3}.prt.model.sel_ops.data_op_mc = D.meancentre;
    if ~isempty(D.data_op)
        if isfield(matlabbatch{3}.prt.model.sel_ops.use_other_ops,'no_op')
            matlabbatch{3}.prt.model.sel_ops.use_other_ops = rmfield(matlabbatch{3}.prt.model.sel_ops.use_other_ops,'no_op');
        end
        matlabbatch{3}.prt.model.sel_ops.use_other_ops.data_op  = D.data_op;
    end
    
    % permutations
    matlabbatch{4}.prt.cv_model.infile = {fullfile(D.spm_path,'PRT.mat')};
    matlabbatch{4}.prt.cv_model.model_name = matlabbatch{3}.prt.model.model_name;
    matlabbatch{4}.prt.cv_model.perm_test=[];
    if D.permtest
        matlabbatch{4}.prt.cv_model.perm_test.perm_t.N_perm = D.permtest;
        matlabbatch{4}.prt.cv_model.perm_test.perm_t.flag_sw = D.saveallweights;
    else
        matlabbatch{4}.prt.cv_model.perm_test.no_perm  = 1;
    end
    
    % weights
    matlabbatch{5}.prt.weights.infile = {fullfile(D.spm_path,'PRT.mat')};
    matlabbatch{5}.prt.weights.model_name = matlabbatch{3}.prt.model.model_name;
    matlabbatch{5}.prt.weights.img_name = '';
    matlabbatch{5}.prt.weights.build_wpr=[];
    if ~isempty(ROIfile)
        matlabbatch{5}.prt.weights.build_wpr.atl_name{1,1} = ROIfile;  
    else
        matlabbatch{5}.prt.weights.build_wpr.no_atl = 0;
    end
    matlabbatch{5}.prt.weights.flag_cwi = D.saveallweights;
end

%if D.para
    save(fullfile(D.spm_path,'matlabbatch'),'matlabbatch');
%elseif D.pronto
%    save(fullfile(D.spm_path,'PRT'),'matlabbatch');
%end
save(fullfile(D.spm_path,'sub_info'),'subID','SubInd');
cd(D.spm_path)
if D.pronto
    if strcmp(D.anatype,'singlesubject')
        Nsub = length(matlabbatch{1}.prt.data.group.select.subject);
        for n = 1:Nsub
            run(n).matlabbatch = matlabbatch;
            run(n).matlabbatch{1}.prt.data.group.select.subject = run(n).matlabbatch{1}.prt.data.group.select.subject(n);
            for c = 1:length(run(n).matlabbatch{3}.prt.model.model_type.classification.class)
                run(n).matlabbatch{3}.prt.model.model_type.classification.class(c).group.subj_nums = 1;
            end
            subpath = fullfile(D.spm_path,['sub' num2str(n)]);
            if ~exist(subpath,'dir')
                mkdir(subpath);
            end
            run(n).matlabbatch{1}.prt.data.dir_name = {subpath};
            run(n).matlabbatch{2}.prt.fs.infile = {fullfile(subpath,'PRT.mat')};
            run(n).matlabbatch{3}.prt.model.infile = {fullfile(subpath,'PRT.mat')};
            run(n).matlabbatch{4}.prt.cv_model.infile = {fullfile(subpath,'PRT.mat')};
            run(n).matlabbatch{5}.prt.weights.infile = {fullfile(subpath,'PRT.mat')};
        end
    else
        run(1).matlabbatch = matlabbatch;
    end
    %jobs = fullfile(D.spm_path,'PRT.mat');
    %inputs = cell(0, 1);
    %job_id = cfg_util('initjob', jobs);
    %sts    = cfg_util('filljob', job_id, matlabbatch{:});
    %if sts
    for n = 1:length(run)
        if n==1
            cfg_get_defaults('cfg_util.genscript_run', @genscript_run);
            cfg_util('initcfg');
            prt_batch
        end
        jid = cfg_util('initjob', run(n).matlabbatch)
        cfg_util('run', jid);
        cfg_util('deljob', jid);
        close all
    end
    %end
    %cfg_util('deljob', job_id);
end
try
    if D.para
        spm_jobman('initcfg')
        spm_jobman('run',matlabbatch);
    end
catch
    if D.para==1
        % reduced number of subjects if fails
        subind = randperm(length(matlabbatch{1, 1}.spm.tools.snpm.des.PairT.fsubject));
        matlabbatch{1, 1}.spm.tools.snpm.des.PairT.fsubject = matlabbatch{1, 1}.spm.tools.snpm.des.PairT.fsubject(subind(1:D.nSubs));
        spm_jobman('initcfg');
        spm_jobman('run',matlabbatch);
    end
end