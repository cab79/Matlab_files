close all
clear all
dbstop if error
% directory containing cosmo projections
cosdir = 'C:\Data\CORE\cosmo';
%cosdir='C:\Data\CORE\ERPs';
% range of cosdir extensions, each with different MM projections
%cosdirext = {'LDA_part4_timechan_0_50','LDA_part4_timechan_50_100','LDA_part4_timechan_100_150','LDA_part4_timechan_150_300'};
cosdirext = {'LDA_part4_timechan_0_150'};
%cosdirext = {'run1_timewin_33_73_basewin_-50_0_centrechan_E93_nNeighbours_2_maskthresh_0.2'};
% generic filename
cosname = 'CORE*_4_mm_projection.mat';
%cosname = '*_mm_proj.mat';

% directory containing HGF projections
hgfdir = 'C:\Data\CORE\Behaviour\July2017\HGF_Results';
%hgfmodels = {'Sim_SDT_RT_a2_noprior','Sim_SDT_RT_a2','Sim_KF_RT_a2','Sim_3lev_RT_a2'};
hgfmodels = {'Sim_SDT_soft_a2_noprior','Sim_SDT_soft_a2','Sim_KF_soft_a2','Sim_3lev_soft_a2'};%,'3lev_Bayes_part4','KF_Bayes_part4'};
%hgfmodels = {'3lev_Bayes_part4'};
%hgfmodels = {'Sim_3lev_RTsoft'};
%hgfmodels = {'Sim_KF-RTsoft'};
hgfname = '_aff.mat';
combine_inputs = 1;

% HGF trajectories
%trajnames = {'mu','sa','muhat','sahat','dau','da','ud','wt','psi','epsi'};
%trajnames = {'dau','da','epsi','ud','null'};
trajnames = {'null1','dau','da'};

results = struct;
for c = 1:length(cosdirext)

    cfiles = dir(fullfile(cosdir,cosdirext{c},cosname));

    for f = 1:length(cfiles)
        cfname = cfiles(f).name;
        load(fullfile(cosdir,cosdirext{c},cfname));
        cosextname = matlab.lang.makeValidName(cosdirext{c});

        % get subject ID and hgf filename for each subject
        Cs = strsplit(cfname,'_');
        sub = strsplit(Cs{1},'CORE');
        sub = sub{2};
        hfname1 = [sub '_bopars' hgfname];
        hfname2 = [sub '_sim' hgfname];

        for m = 1:length(hgfmodels)
            try
                mod = load(fullfile(hgfdir,hgfmodels{m},hfname1));
            catch
                try
                    mod = load(fullfile(hgfdir,hgfmodels{m},hfname2));
                catch
                    continue
                end
            end

            % make valid fieldname
            modname = hgfmodels{m};
            modname = matlab.lang.makeValidName(modname);

            % in case the model has bopars instead of sim, make it sim
            mfield = fieldnames(mod);
            sim = mod.(mfield{1});

            % find inputs
            in = sim.u(:,1);
            incond = sim.u(:,1) .* sim.u(:,2);
            if combine_inputs
                ni = 1;
                iii{1}=tnums;
                IB{1} = 1:length(tnums);
            else
                n_in = unique(in);
                n_in = n_in(~isnan(n_in));
                ni = length(n_in);
                for i = 1:ni
                    ii{i} = find(in==n_in(i));
                    [iii{i},~,IB{i}] = intersect(ii{i},tnums);
                end
            end

            % for each input type
            for i = 1:ni
                try
                    mmi = zscore(mm(IB{i}));
                catch
                    mmi = zscore(MM.mm(IB{i}))';
                end

                trcount=0;
                for tr = 1:length(trajnames)
                    trajname = trajnames{tr};
                    if strcmp(trajname,'null1')
                        trajmat = in;
                    else
                        trajmat = sim.traj.(trajname);
                    end
                    for tl = 1:size(trajmat,2)
                        if all(isnan(trajmat(:,tl)))
                            continue
                        end
                        trcount = trcount+1;
                        traj = trajmat(:,tl);
                        traj = zscore(traj(iii{i}));

                        [r,p] = corr(traj,mmi','type','Spearman');
                        %[r,p] = corr(traj,mmi','type','Pearson');

                        results.(cosextname).(modname).([trajname num2str(tl)]).r(f,i) = r;
                        results.(cosextname).(modname).([trajname num2str(tl)]).p(f,i) = p;

                        %% spm_peb
                        if ~isnan(r)
                            clear C P F
                            P{1}.X = traj;
                            %settings of http://www.sciencedirect.com/science/article/pii/S1053811912004466#f0015
                            % Cb is prior for 2nd level parameters and in spm_peb comes from Eq. 16 in this paper and is all zeros
                            % except for the final element of exp(32)
                            %leads to "singular or badly scaled matrix" to include level 2, and
                            %results are the same as without 2nd level
                            %P{2}.X = 0;
                            %P{2}.C{1} =1;
                            [C,P,F] = spm_PEB(mmi',P);
                            results.(cosextname).(modname).([trajname num2str(tl)]).F(f,i) = F;
                        else
                            results.(cosextname).(modname).([trajname num2str(tl)]).F(f,i) = -10000;
                        end
                        %%
                        %figure
                        %scatter(traj,mmi')
                        %lsline
                        %figure
                        %scatter(tiedrank(traj),tiedrank(mmi'))
                        %lsline
                    end
                end
            end
        end
    end
end
savename = ['Traj_corr_results_' datestr(now,30)];
save(fullfile(cosdir,savename),'results','ni');

% combine results
if 0
    file1 = 'C:\Data\CORE\cosmo\LDA_part4_timechan\Null_model_results.mat';
    file2 = ['C:\Data\CORE\cosmo\LDA_part4_timechan\' savename];
    f1=load(file1);
    f2=load(file2);
    fieldNames = fieldnames(f1.results);
    for i = 1:size(fieldNames,1)
        f2.results.(fieldNames{i}) = f1.results.(fieldNames{i});
    end
    results = f2.results;
end

%% Bayesian Model Comparison
% for each input type, which HGF trajectory is most closely related to the EEG projection over all subjects?

close all
clear Fs
%select cosmo directory or leave blank to analyse all in the results file
cosmoselect = '';
%model comparison for a single trajectory
modelselect = ''; % blank to compare all models
%trajectory comparison for a single model
trajselect = ''; % blank to compare all trajectories

% restructure data for VBA
cosmonames = fieldnames(results);
if ~isempty(cosmoselect); cosmoselect = matlab.lang.makeValidName(cosmoselect);end;
nmod = 0;
names = {};
for c = 1:length(cosmonames)
    % select cosmo trajectory
    if ~isempty(cosmoselect)
        if ~strcmp(cosmoselect,cosmonames{c})
            continue
        end
    end
    modelnames = fieldnames(results.(cosmonames{c}));
    if ~isempty(modelselect); modelselect = matlab.lang.makeValidName(modelselect);end;
    
    for m = 1:length(modelnames)
        % select model
        if ~isempty(modelselect)
            if ~strcmp(modelselect,modelnames{m})
                continue
            end
        end
        trajnames = fieldnames(results.(cosmonames{c}).(modelnames{m}));
        for tr = 1:length(trajnames)

            % select trajectory
            if ~isempty(trajselect)
                if ~strcmp(trajselect,trajnames{tr})
                    continue
                end
            end
            nmod = nmod+1;

            % extract data
            Ft=results.(cosmonames{c}).(modelnames{m}).(trajnames{tr}).F;
            for nin = 1:ni; % for each input type
                Fs{nin}(:,nmod) = Ft(:,nin);
                names{length(names)+1,1} = [cosmonames{c} '_' modelnames{m} '_' trajnames{tr} '_' num2str(nin)];
            end
        end
    end
end

% remove rows with zero values
rem_rows = find(any((Fs{nin}==0)'));
Fs{nin}(rem_rows,:)=[];

% remove columns with low evidence
rem_cols = find(any(Fs{nin}==-10000));
Fs{nin}(:,rem_cols)=[];

% apply VBA
for nin = 1:ni; % for each input type
    [VBAposterior,VBAout] = VBA_groupBMC(Fs{nin}')
end
names


compare=[1 2];
% Bayes factors
BF = exp(Fs{nin}(:,compare(1)) - Fs{nin}(:,compare(2)));
% Group Bayes factor - sensitive to outliers!
GBF = prod(BF);
