function S=erp_freq_analysis(S)
%% ERP, Frequency, Time-Frequency and Coherence analysis
% To be run after data can be cleaned 
% Before running this, check data quality in EEGLAB
% Requires Fieldtrip to be in the Matlab path

S.nMarkType = length(S.markers);
if strcmp(S.analysistype,'ERP')
%    dp = 250;
    S.nf = 1;
%    times = -0.2:0.004:0.796;
elseif strcmp(S.analysistype,'TF')
%    dp = 181;
    S.nf = length(S.freqsrange);
%    intimes = -0.2:0.004:0.796;
end
if strcmp(S.analysistype,'ERP') && S.CSD.apply 
    load(S.CSD.montage); % path to CSD montage (only needed if using CSD)
end 

% GET FILE LIST
S.filepath = S.setpath;
S = getfilelist(S);

for f = 1:length(S.filelist)
    file = S.filelist{f};
    EEG = pop_loadset('filename',file,'filepath',S.setpath);

    % SAVE CHANLOCS FOR PLOTTING LATER (needed for EEGLAB's topoplot function)
    if f==1
        if ~exist(fullfile(S.setpath,'chanlocs.mat'),'file')
            chanlocs = EEG.chanlocs;
            save(fullfile(S.setpath,'chanlocs.mat'),'chanlocs')
        end
    end

    % Remove ECG if not done already
    ecgchan = strcmp({EEG.chanlocs.labels},'ECG');
    if any(ecgchan)
        EEG = pop_select(EEG,'nochannel',find(ecgchan));
    end
    
    % remove baseline
    if S.rmbase
        EEG = pop_rmbase( EEG, [S.basewin(1)*1000 S.basewin(2)*1000]);
    end
    
    dsize = size(EEG.data);
    S.intimes=EEG.times;
    %SINGDAT = cell(1);
    
    % prepare for saving
    [pth nme ext] = fileparts(file); 
    sname_ext = [S.analysistype '.mat'];
    sname = [nme '_' sname_ext];
    
    % switch to selected analysis
    switch S.analysistype
        case 'ERP'
            EEGall=EEG;
            for mt = 1:S.nMarkType
                if S.combinemarkers
                    EEG = pop_selectevent(EEGall,'type',S.markers);
                else
                    EEG = pop_selectevent(EEGall,'type',S.markers{mt});
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
                if S.CSD.apply
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
            if ~exist(S.erppath,'dir')
                mkdir(S.erppath);
            end
            save(fullfile(S.erppath,sname),'tldata');
            
        case {'Freq','TF','Coh'}
            
            % convert to Fieldtrip
            FTEEG=convertoft(EEG);
            FTEEG.event=EEG.event;
            
            if strcmp(S.freqtype,'induced') 
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
            if S.mov_avg_trials
                nAvg = floor((size(fdata.powspctrm,1)/S.mov_avg_trials_step)-(S.mov_avg_trials/S.mov_avg_trials_step))+1;
                for ma = 1:nAvg
                    st_ind = (ma-1)*S.mov_avg_trials_step+1;
                    fdata.madata(:,:,ma) = squeeze(nanmean(fdata.powspctrm(st_ind:st_ind+S.mov_avg_trials-1,:,:),1));
                end
            end

            % COHERENCE (WPLI)
            if strcmp(S.analysistype,'Coh')
                % coherence matrices
                nchan = length(fdata.label);
                matrix=zeros(size(S.freqsrange,1),nchan,nchan); 
                coh = zeros(nchan,nchan);
                cohboot = zeros(nchan,nchan);

                wpli = ft_connectivity_wpli(fdata.crsspctrm,'debias',true,'dojack',false);
                for c = 1:size(S.freqsrange,1)
                    fprintf('f %d',c);
                    [M, bstart] = min(abs(fdata.freq-freqlist(c,1)));
                    [M, bend] = min(abs(fdata.freq-freqlist(c,2)));
                    coh(:) = 0;
                    coh(logical(tril(ones(size(coh)),-1))) = max(wpli(:,bstart:bend),[],2);
                    coh = tril(coh,1)+tril(coh,1)';
                    matrix(c,:,:) = coh;
                end

                for nboot = 1:S.bootrep
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
            freqs = S.freqsrange;

            if strcmp(S.analysistype,'TF') && ~isempty(S.baselinetype)
                % remove baseline
                cfg.baseline     = S.basewin;
                cfg.baselinetype = S.baselinetype; %'absolute', 'relative', 'relchange', 'normchange' or 'db' (default = 'absolute')
                cfg.parameter    = 'powspctrm';
                fdata = ft_freqbaseline_CAB(cfg, fdata);
            end
            
            % separate marker types
            if ~S.combinemarkers
                events = {FTEEG.event.type};
                for mt = 1:S.nMarkType
                    cfg.trials      = find(strcmp(events,S.markers{mt}));
                    cfg.avgoverrpt  = 'no';
                    fdatat{mt} = ft_selectdata(cfg, fdata);
                end
                fdata=fdatat;
                clear fdatat;
            end
            
            % SAVE
            if ~exist(S.freqpath,'dir')
                mkdir(S.freqpath)
            end
            save(fullfile(S.freqpath,sname),'fdata');
    end
    
end