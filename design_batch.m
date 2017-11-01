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
%%
%-----------------------------------------------------------------------  
function D=design_batch(D)
dbstop if error

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


% select parametric or non-parametric analyses and the range of models
if D.para==1
    load(D.ffbatch);
    paraname = '_spm';
elseif D.para==2
    load(D.npbatch);
    paraname = '_snpm';
    ismain = any(D.maineffects);
    isinter = any(D.interactions);
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

factnames = '';
for f = 1:length(D.factors)
    factnames = [factnames '_' D.factors{f}];
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

[~,~,pdata] = xlsread(D.pdatfile);
grp_col = find(strcmp(pdata(1,:),D.grphead));
sub_col = find(strcmp(pdata(1,:),D.subhead));
inc_col = find(strcmp(pdata(1,:),D.inchead));

if isempty(grp_col) && ~isfield(D,'grpind')
    error(['no header in xls file called ' D.grphead]);
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
    %Change char Grp inputs to numbers
    grpdat = pdata(2:end,grp_col);
    if isnumeric(grpdat{2,1})
        grptype = unique([grpdat{:}]);
        grptype = grptype(~isnan(grptype));
        grptype(grptype==0)=[];
        Ngrp = length(grptype);
    else
        grptype = unique(grpdat);
        grptype(isempty(grptype))=[];
        Ngrp = length(grptype);
        for g = 1:Ngrp
            grp_idx = cellfun(@(x) any(strcmp(grptype(g),x)), grpdat, 'UniformOutput', 0);
            grpdat(cell2mat(grp_idx)) = {[g]};
        end
    end
end

% find index of subjects to include in analysis
SubInd = cell(Ngrp,1);
Subs = [];
inc_idx = cellfun(@(x) ismember(x,D.include_codes), pdata(2:end,inc_col), 'UniformOutput', 0);
inc_idx = find(cell2mat(inc_idx));

% find subject indices for each specific group
for g = 1:Ngrp
    grp_idx = find(cellfun(@(x) x==g, grpdat, 'UniformOutput', 1));
    if D.useIDfile==1
        SubInd{g,1} = intersect(inc_idx,grp_idx);
    else
        SubInd{g,1} = grp_idx;
    end
    Nsub(g,1) = length(SubInd{g,1});
end

% add factors to the design
if D.para==1
    matlabbatch{1}.spm.stats.factorial_design.dir = {D.spm_path};
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
elseif D.para==2
    generic.dir = {D.spm_path};
    generic.nPerm = D.nPerm;
    generic.vFWHM = D.vFWHM;
    generic.bVolm = D.bVolm;
    generic.ST.ST_U = D.ST_U;
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

gs = 0;
subID={};
for g = 1:Ngrp
    for s = 1:Nsub(g)
        gs = gs+1;
        if D.useIDfile==1
            subID{gs} = deblank(pdata{SubInd{g,1}(s)+1,sub_col});
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
        
        subimg = cell(1,1);
        if D.para==1
            % for spm, load all files
            for i = 1:length(D.imglist)
                fname = fullfile(subdir, [subfile D.imglist{i}]);
                if ~exist(fname,'file')
                    fname_temp = dir(fullfile(subdir, [subfile D.imglist{i}]));
                    if length(fname_temp)>1
                        display(fname);
                        error('Subject file name is not unique in this folder');
                    elseif length(fname_temp)==0
                        display(fname);
                        error('no file found with this name');
                    end
                    fname=fullfile(subdir, fname_temp.name);
                end
                subimg{i,1} = [fname ',1'];
            end
        elseif D.para==2
            if D.useIDfile==1
                % for SnPM, check averaged files exist; if not then create them.
                fnames = dir(fullfile(subdir, ['*' scanname '*']));
                if wnum==0
                    condlist = ones(length(D.cond_list),1);
                else
                    condlist = D.cond_list(:,wnum);
                end
                if isempty(fnames)
                    factor_img(subdir,D.imglist,condlist,scanname);
                    fnames = dir(fullfile(subdir, ['*' scanname '*']));
                end
                for i=1:length(fnames)
                    subimg{i,1} = fullfile(subdir,[fnames(i).name ',1']);
                end
            else
                subimg = subfile;
            end
        end
        scans = subimg;
        
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
        
        % add the above scan and group info to the matlabbatch, or create new scans
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
        end
        
    end
end

% add main effects and interactions
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

% covariates
if ~isempty(D.cov_names)
    if ~strcmp(D.cov_names{1},'')
        for c = 1:length(D.cov_names)
            cov_col = find(strcmp(pdata(1,:),D.cov_names{c}));
            covdat = pdata(2:end,cov_col);

            % convert categorical (char) covariate data to numbers
            if ~isnumeric(covdat{2,1})
                covtype = unique(covdat);
                Ncovtype = length(covtype);
                for c = 1:Ncovtype
                    cov_idx = cellfun(@(x) any(strcmp(covtype(c),x)), covdat, 'UniformOutput', 0);
                    covdat(cell2mat(cov_idx)) = {[c]};
                end
            end

            covdat = [covdat{:}]';
            covdat = covdat(vertcat(SubInd{:}));
            covdat = repmat(covdat(inc_idx),1,length(scans));
            covdat = reshape(covdat', size(covdat,1)*size(covdat,2),1);
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
        end
    end
else
    if D.para==1
        matlabbatch{1}.spm.stats.factorial_design.cov = struct('c', {}, 'cname', {}, 'iCFI', {}, 'iCC', {}); 
    elseif D.para==2
        if strcmp(model,'g')
            matlabbatch{1}.spm.tools.snpm.des.TwoSampT.cov = struct('c', {}, 'cname', {}); 
        elseif strcmp(model,'gw')
            matlabbatch{1}.spm.tools.snpm.des.TwoSampPairT.cov = struct('c', {}, 'cname', {});
        end
    end
end
if D.para==1
    matlabbatch{1}.spm.stats.factorial_design.multi_cov = struct('files', {}, 'iCFI', {}, 'iCC', {});
end

% specify mask
if ~isempty(D.time_ana)
    maskname = [num2str(D.time_ana(1)) '_' num2str(D.time_ana(2)) '.img'];
    maskfile = fullfile(D.mask_path,maskname);
    if ~exist(maskfile,'file')
        S = struct;
        S.image = subimg{1}(1:end-2); % example image
        S.timewin = D.time_ana;
        S.outfile = maskfile;
        spm_eeg_mask(S)
    end
    masking.em = {[maskfile  ',1']};
else
    masking.em = {''};
end
masking.tm.tm_none = 1;
masking.im = 1;
if D.para==1
    matlabbatch{1}.spm.stats.factorial_design.masking = masking;
    matlabbatch{1}.spm.stats.factorial_design.globalc.g_omit = 1;
    %globalm=matlabbatch{1}.spm.stats.factorial_design.globalm;
elseif D.para==2
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
end


%specify scaling parameters
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
    if strcmp(model,'w') || strcmp(model,'ww')
        matlabbatch{1}.spm.tools.snpm.des.PairT.globalm = globalm;
    elseif strcmp(model,'g')
        matlabbatch{1}.spm.tools.snpm.des.TwoSampT.globalm = globalm;
    elseif strcmp(model,'gw')
        matlabbatch{1}.spm.tools.snpm.des.TwoSampPairT.globalm = globalm;
    elseif strcmp(model,'one')
        matlabbatch{1}.spm.tools.snpm.des.OneSampT.globalm = globalm;
    end

    matlabbatch{2}.spm.tools.snpm.cp.snpmcfg = {fullfile(D.spm_path,'SnPMcfg.mat')};

    matlabbatch{3} = struct();
    matlabbatch{3}.spm.tools.snpm.inference.SnPMmat = {fullfile(D.spm_path,'SnPM.mat')};
    matlabbatch{3}.spm.tools.snpm.inference.Thr.Clus.ClusSize.CFth = NaN; % set to NaN if already defined in the design
    matlabbatch{3}.spm.tools.snpm.inference.Thr.Clus.ClusSize.ClusSig.FWEthC = 0.05;
    matlabbatch{3}.spm.tools.snpm.inference.Tsign=1;
    matlabbatch{3}.spm.tools.snpm.inference.WriteFiltImg.name='SnPM_filtered';
    matlabbatch{3}.spm.tools.snpm.inference.Report = 'MIPtable';
end

save(fullfile(D.spm_path,'matlabbatch'),'matlabbatch');
save(fullfile(D.spm_path,'sub_info'),'subID','SubInd');
spm_jobman('initcfg')
try
    spm_jobman('run',matlabbatch);
catch
    % reduced number of subjects if fails
    subind = randperm(length(matlabbatch{1, 1}.spm.tools.snpm.des.PairT.fsubject));
    matlabbatch{1, 1}.spm.tools.snpm.des.PairT.fsubject = matlabbatch{1, 1}.spm.tools.snpm.des.PairT.fsubject(subind(1:D.nSubs));
    spm_jobman('initcfg');
    spm_jobman('run',matlabbatch);
end