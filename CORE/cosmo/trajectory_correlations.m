close all
clear all
dbstop if error
% directory containing cosmo projections
cosdir = 'C:\Data\CORE\cosmo';
% range of cosdir extensions, each with different MM projections
%cosdirext = {'LDA_part4_timechan_0_50','LDA_part4_timechan_50_100','LDA_part4_timechan_100_150','LDA_part4_timechan_150_300'};
cosdirext = {'LDA_part4_timechan'};
% generic filename
cosname = 'CORE*_4_mm_projection.mat';

% directory containing HGF projections
hgfdir = 'C:\Data\CORE\Behaviour\July2017\HGF_Results';
%hgfmodels = {'Sim_KF-RT','Sim_KF-RTsoft','Sim_KF-soft','Sim_3lev-RT','Sim_3lev_RTsoft','Sim_3lev-soft','3lev_Bayes_part4','KF_Bayes_part4'};
hgfmodels = {'Sim_3lev_RTsoft','Sim_KF-RT','3lev_Bayes_part4','KF_Bayes_part4'};
%hgfmodels = {'Sim_3lev_RTsoft'};
%hgfmodels = {'Sim_KF-RTsoft'};
hgfname = '_aff.mat';
combine_inputs = 1;

% HGF trajectories
%trajnames = {'mu','sa','muhat','sahat','dau','da','ud','wt','psi','epsi'};
trajnames = {'dau','da','epsi','ud'};

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
                mod = load(fullfile(hgfdir,hgfmodels{m},hfname2));
            end

            % make valid fieldname
            modname = hgfmodels{m};
            modname = matlab.lang.makeValidName(modname);

            % in case the model has bopars instead of sim, make it sim
            mfield = fieldnames(mod);
            sim = mod.(mfield{1});

            % find inputs
            in = sim.u(:,1);
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
                mmi = zscore(mm(IB{i}));

                trcount=0;
                for tr = 1:length(trajnames)
                    trajname = trajnames{tr};
                    trajmat = sim.traj.(trajname);
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
                            results.(cosextname).(modname).([trajname num2str(tl)]).F(f,i) = -100000;
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
save(fullfile(cosdir,['Traj_corr_results_' datestr(now,30)]),'results','ni');


% combine results
if 0
    file1 = 'C:\Data\CORE\cosmo\LDA_part4_timechan\Traj_corr_results_20170729T203536.mat';
    file2 = 'C:\Data\CORE\cosmo\LDA_part4_timechan\Traj_corr_results_20170729T195853.mat';
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
trajselect = 'da1'; % blank to compare all trajectories

% restructure data for VBA
cosmonames = fieldnames(results);
if ~isempty(cosmoselect); cosmoselect = matlab.lang.makeValidName(cosmoselect);end;
nmod = 0;
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
            for in = 1:ni; % for each input type
                Fs{in}(:,nmod) = Ft(:,in);
            end
        end
    end
end

% apply VBA
for in = 1:ni; % for each input type
    [VBAposterior,VBAout] = VBA_groupBMC(Fs{in}')
end
