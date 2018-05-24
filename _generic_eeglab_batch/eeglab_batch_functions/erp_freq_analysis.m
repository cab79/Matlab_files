function S=erp_freq_analysis(S)
%% ERP, Frequency, Time-Frequency and Coherence analysis
% To be run after data can be cleaned 
% Before running this, check data quality in EEGLAB
% Requires Fieldtrip to be in the Matlab path
S.func = 'tf';
S.(S.func).nMarkType = length(S.(S.func).epoch.markers);
if strcmp(S.(S.func).select.datatype,'ERP')
%    dp = 250;
    S.(S.func).nf = 1;
%    times = -0.2:0.004:0.796;
elseif strcmp(S.(S.func).select.datatype,'TF')
%    dp = 181;
    S.(S.func).nf = length(S.(S.func).select.freq);
%    intimes = -0.2:0.004:0.796;
end
if strcmp(S.(S.func).select.datatype,'ERP') && S.(S.func).CSD.apply 
    load(S.(S.func).CSD.montage); % path to CSD montage (only needed if using CSD)
end 

% GET FILE LIST
S.path.file = S.path.prep;
S = getfilelist(S);

for f = 1:length(S.(S.func).filelist)
    file = S.(S.func).filelist{f};
    EEG = pop_loadset('filename',file,'filepath',S.path.prep);

    % SAVE CHANLOCS FOR PLOTTING LATER (needed for EEGLAB's topoplot function)
    %if f==1
    %    if ~exist(fullfile(S.path.prep,'chanlocs.mat'),'file')
    %        chanlocs = EEG.chanlocs;
    %        save(fullfile(S.path.prep,'chanlocs.mat'),'chanlocs')
    %    end
    %end
    
    % add in any missing chans, fill with nan
    load(fullfile(S.path.prep,'chanlocs.mat'));
    labs = {chanlocs.labels};
    actlabs = {EEG.chanlocs.labels};
    missing = find(~ismember(labs,actlabs));
    if ~isempty(missing)
        for m = 1:length(missing);
            EEG.chanlocs(missing(m)+1:end+1) = EEG.chanlocs(missing(m):end);
            EEG.chanlocs(missing(m)) = chanlocs(missing(m));
            EEG.data(missing(m)+1:end+1,:,:) = EEG.data(missing(m):end,:,:);
            EEG.data(missing(m),:,:) = nan;
        end
        EEG.nbchan = length(EEG.chanlocs);
        if EEG.nbchan~=length(labs)
            error('not enough chans')
        end
    end

    % Remove ECG if not done already
    ecgchan = strcmp({EEG.chanlocs.labels},'ECG');
    if any(ecgchan)
        EEG = pop_select(EEG,'nochannel',find(ecgchan));
    end
    
    % remove baseline
    if S.(S.func).epoch.rmbase
        EEG = pop_rmbase( EEG, [S.(S.func).epoch.basewin(1)*1000 S.(S.func).epoch.basewin(2)*1000]);
    end
    
    dsize = size(EEG.data);
    S.(S.func).intimes=EEG.times;
    %SINGDAT = cell(1);
    
    % prepare for saving
    [pth nme ext] = fileparts(file); 
    sname_ext = [S.(S.func).select.datatype '.mat'];
    sname = [nme '_' sname_ext];
    
    % switch to selected analysis
    switch S.(S.func).select.datatype
        case 'ERP'
            EEGall=EEG;
            for mt = 1:S.(S.func).nMarkType
                if S.(S.func).epoch.combinemarkers
                    EEG = pop_selectevent(EEGall,'type',S.(S.func).epoch.markers);
                else
                    try
                        EEG = pop_selectevent(EEGall,'type',S.(S.func).epoch.markers{mt});
                    catch
                        continue
                    end
                end
                
                % convert to Fieldtrip
                FTEEG=convertoft(EEG);
                FTEEG.event=EEG.event;
                
                % ERP: timelocked (evoked) data
                cfg = [];
                cfg.keeptrials = 'yes';
                cfg.feedback = 'textbar';
                tldata{mt} = ft_timelockanalysis(cfg,FTEEG);

                % apply CSD
                if S.(S.func).CSD.apply
                    [ntrial,nchan,nsamp] = size(tldata{mt}.trial);
                    trials = reshape(permute(tldata{mt}.trial,[2 3 1]),nchan,nsamp*ntrial);
                    trials=CSD(trials,G,H); 
                    trials = permute(reshape(trials,nchan,nsamp,ntrial),[3 1 2]);
                    tldata{mt}.trial = trials;
                end
                tldata{mt} = ft_timelockanalysis(cfg,tldata{mt});
                %erp = squeeze(double(mean(tldata.trial,1)));
            end
            
            % SAVE
            if ~exist(S.path.erp,'dir')
                mkdir(S.path.erp);
            end
            save(fullfile(S.path.erp,sname),'tldata');
            
        case {'Freq','TF','Coh'}
            
            % convert to Fieldtrip
            FTEEG=convertoft(EEG);
            FTEEG.event=EEG.event;
            
            if strcmp(S.(S.func).freq.type,'induced') 
                % timelocked (evoked) data
                cfg = [];
                cfg.keeptrials = 'yes';
                cfg.feedback = 'textbar';
                tldata = ft_timelockanalysis(cfg,FTEEG);
                % subtract out ERP for induced activity
                for t = 1:length(FTEEG)
                    FTEEG.trial{t} = FTEEG.trial{t}-tldata.avg; 
                end
            end
            
            % freq analysis
            fdata = FTfreqanalysis(FTEEG,S);
            fdata.time = FTEEG.time;
            act_freq = fdata.freq;
            disp(['actual frequencies: ' num2str(act_freq)])
            
            % moving averages over time
            if S.(S.func).op.mov_avg_trials
                nAvg = floor((size(fdata.powspctrm,1)/S.(S.func).op.mov_avg_trials_step)-(S.(S.func).op.mov_avg_trials/S.(S.func).op.mov_avg_trials_step))+1;
                for ma = 1:nAvg
                    st_ind = (ma-1)*S.(S.func).op.mov_avg_trials_step+1;
                    fdata.madata(:,:,ma) = squeeze(nanmean(fdata.powspctrm(st_ind:st_ind+S.(S.func).op.mov_avg_trials-1,:,:),1));
                end
            end

            % COHERENCE (WPLI)
            if strcmp(S.(S.func).select.datatype,'Coh')
                % coherence matrices
                nchan = length(fdata.label);
                matrix=zeros(size(S.(S.func).select.freq,1),nchan,nchan); 
                coh = zeros(nchan,nchan);
                cohboot = zeros(nchan,nchan);

                wpli = ft_connectivity_wpli(fdata.crsspctrm,'debias',true,'dojack',false);
                for c = 1:size(S.(S.func).select.freq,1)
                    fprintf('f %d',c);
                    [M, bstart] = min(abs(fdata.freq-freqlist(c,1)));
                    [M, bend] = min(abs(fdata.freq-freqlist(c,2)));
                    coh(:) = 0;
                    coh(logical(tril(ones(size(coh)),-1))) = max(wpli(:,bstart:bend),[],2);
                    coh = tril(coh,1)+tril(coh,1)';
                    matrix(c,:,:) = coh;
                end

                for nboot = 1:S.(S.func).freq.bootrep
                    fprintf('nboot %d',nboot);
                    for tr = 1:length(fdata.trial)
                        for ele = 1:size(fdata.trial{tr},1);
                            trial = fdata.trial{tr};
                            surrEEG=[];
                            surrEEG = [phaseran(squeeze(trial(tr,:)'),1); 0]';
                            %surrEEG = [phaseran(EEG.trial{1}',1); zeros(1, length(EEG.label))]';
                            trial(tr,:) = surrEEG;
                        end
                        fdata.trial{tr} = trial;
                    end

                    fdata = FTfreqanalysis(fdata,S);
                    wpli_boot = ft_connectivity_wpli(fdata.crsspctrm,'debias',true,'dojack',false);

                    for c = 1:size(freqlist,1)
                        [M, bstart] = min(abs(fdata.freq-freqlist(c,1)));
                        [M, bend] = min(abs(fdata.freq-freqlist(c,2)));

                        cohboot(:) = 0;
                        cohboot(logical(tril(ones(size(cohboot)),-1))) = max(wpli_boot(:,bstart:bend),[],2);
                        cohboot = tril(cohboot,1)+tril(cohboot,1)';
                        bootmat(c,:,:,nboot) = cohboot;
                    end
                    clear wpli_boot
                end
            end
            freqs = S.(S.func).select.freq;

            if strcmp(S.(S.func).select.datatype,'TF') && ~isempty(S.(S.func).freq.basenorm)
                % remove baseline
                cfg.baseline     = S.(S.func).epoch.basewin;
                cfg.baselinetype = S.(S.func).freq.basenorm; %'absolute', 'relative', 'relchange', 'normchange' or 'db' (default = 'absolute')
                cfg.parameter    = 'powspctrm';
                fdata = ft_freqbaseline_CAB(cfg, fdata);
            end
            
            % separate marker types
            if ~S.(S.func).epoch.combinemarkers
                % get marker occuring at zero latency
                events={};
                for m=1:length(EEG.epoch)
                    events{m} = EEG.epoch(m).eventtype{find([EEG.epoch(m).eventlatency{:}]==0)};
                end
                for mt = 1:S.(S.func).nMarkType
                    cfg.trials      = find(strcmp(events,S.(S.func).epoch.markers{mt}));
                    cfg.avgoverrpt  = 'no';
                    fdatat{mt} = ft_selectdata(cfg, fdata);
                end
            else
                fdatat = {fdata};
            end
            
            if isfield(fdata,'madata')
                fdatat{1}.madata=fdata.madata;
            end
            
            fdata=fdatat;
            clear fdatat;
            % SAVE
            switch S.(S.func).select.datatype
                case 'Freq'
                    if ~exist(S.path.freq,'dir')
                        mkdir(S.path.freq)
                    end
                    save(fullfile(S.path.freq,sname),'fdata');
                case 'TF'
                    if ~exist(S.path.tf,'dir')
                        mkdir(S.path.tf)
                    end
                    save(fullfile(S.path.tf,sname),'fdata');
            end
                    
    end
    
end