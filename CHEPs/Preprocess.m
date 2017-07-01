%% PREPROCESSING
% 1. filtering, epoching, baseline correction, rereference
% 2. rejects chans/trials only with high frequency noise
% 3. reject ICA components related to eye movement (and 50Hz if no notch applied) but not that related to individual channel noise 
% 4. reject remaining chans/trials outside normal range
%%

clear all
origpath = 'C:\Data\CHEPs\Orig';
anapath = 'C:\Data\CHEPs\Preprocessed';
cd(origpath);
files = dir('*orig.set');
cd(anapath);
combine_all=0; % combining left and right stimulations, or that of different experiments, may be unwise for ICA purposes.

timebins = [-0.5 1]; % for epoching, TSOT(2)
filterset = [0.5 40];
notch_on = 1;
addpath(genpath('M:\Matlab\Matlab_files\CORE\Supporting functions'));

files_ana = [1]% 1:length(files);
markers = {'STIM' 'DIN2'};
use_marker = 2;
STIM_delay = 0;

for f = files_ana
    
    orig_file = files(f).name;
    [pth nme ext] = fileparts(orig_file); 
    C = strsplit(nme,'_');
    
    %if strcmp(C{2},'2')
        timebin = timebins(1,:);
    %elseif strcmp(C{2},'4')
    %    timebin = timebins(2,:);
    %    isi = ISIs(2);
    %end
    
    EEG = pop_loadset('filename',orig_file,'filepath',origpath);
    
    % FILTER
    if filterset(1)>0; EEG = pop_eegfiltnew( EEG, filterset(1), 0, [], 0);end
    if filterset(2)>0; EEG = pop_eegfiltnew( EEG, 0, filterset(2), [], 0);end
    if notch_on==1 
        EEG = pop_eegfiltnew(EEG,45,55,[],1); 
        EEG = pop_eegfiltnew(EEG,95,105,[],1); 
    elseif notch_on==2
        EEG.data = rm50Hz(EEG.data,[2 1 3],EEG.srate,20);
    end;
    
    if strcmp(markers{use_marker},'STIM')
        for i=1:length(EEG.event)
            if strcmp(EEG.event(i).type,'STIM')
                EEG.event(i).init_time = EEG.event(i).init_time-STIM_delay; 
                EEG.event(i).latency = EEG.event(i).latency-round(STIM_delay*EEG.srate); 
            end
        end
    end
    
    % EPOCH
    EEG = pop_epoch( EEG, markers(use_marker), timebin, 'newname', [C{1} '_' C{2} '_epochs'], 'epochinfo', 'yes');
    
    % Correct STIMs
    %if strcmp(markers{use_marker},'DIN2')
    %    for ep = 1:length(EEG.epoch)
    %        stimevidx = find(strcmp('STIM',EEG.epoch(ep).eventtype));
    %        if ep<length(EEG.epoch); stimevidx1 = find(strcmp('STIM',EEG.epoch(ep+1).eventtype));end;
    %        if isempty(stimevidx)
    %            if length(stimevidx1)>1
    %                event = EEG.epoch(ep+1).eventurevent{1, stimevidx1(1)};
    %                event_idx = find([EEG.event.urevent]==event);
    %                %EEG = pop_editeventvals(EEG, 'changefield', {event_idx 'latency' eeg_point2lat(EEG.event(event_idx-1).latency,[],EEG.srate, [EEG.xmin EEG.xmax])});
    %                EEG.event(event_idx).latency = EEG.event(event_idx-1).latency;
    %            else
    %                error(['not enough stims in epoch ' num2str(ep+1)]);
    %            end
    %        end
    %    end
    %    % EPOCH AGAIN AFTER STIM CORRECTION
    %    EEG = pop_epoch( EEG, {'DIN2'}, timebin, 'newname', [C{1} '_' C{2} '_epochs'], 'epochinfo', 'yes');
    %end


    %EEG = rmTrailsISIvar(EEG,0.003,timebin,isi,markers{use_marker});
        
    % RE-REFERENCE
    EEG = pop_reref( EEG, []);
    
    % LINEAR DETREND
    for i = 1:EEG.trials, EEG.data(:,:,i) = detrend(EEG.data(:,:,i)')'; end;
    
    % REMOVE BASELINE
    EEG = pop_rmbase( EEG, [timebin(1)*100    0]);
    
    
    sname = [C{1} '_' C{2} '_epoched.set'];
    EEG = pop_saveset(EEG,'filename',sname,'filepath',anapath); 
    
end

for f = files_ana
    orig_file = files(f).name;
    [pth nme ext] = fileparts(orig_file); 
    C = strsplit(nme,'_');
    lname = [C{1} '_' C{2} '_epoched.set'];
    EEG = pop_loadset('filename',lname,'filepath',anapath);
    
    % strategy - only remove a very small number of very bad trials / chans
    EEG = FTrejman(EEG,[20 30]); % high freq to identify noise not related to eye movement
    EEG = FTrejman(EEG,[0 5]); % low freq to identify eye movement
    
    %EEG = pop_eegplot(EEG);
    
    sname = [C{1} '_' C{2} '_manrej.set'];
    EEG = pop_saveset(EEG,'filename',sname,'filepath',anapath); 
end

for f = files_ana

    orig_file = files(f).name;
    [pth nme ext] = fileparts(orig_file); 
    C = strsplit(nme,'_');
    lname = [C{1} '_' C{2} '_manrej.set'];
    EEG = pop_loadset('filename',lname,'filepath',anapath);
    
    %[conds, tnums, fnums, bnums] = get_markers(EEG);
    
    %EEGall = EEG;
    %for hand = 1:2
        
        %if hand == 1
        %    markers = leftmarkers;
        %    handname = 'left';
        %else
        %    markers = rightmarkers;
        %    handname = 'right';
        %end
        
        % create an index of conds
        %selectepochs=[];
        %for i = 1:length(markers)
        %    selectepochs = [selectepochs find(conds==markers(i))];
        %end
        %selectepochs = sort(selectepochs);
    
        %EEG = pop_select(EEGall,'trial',selectepochs);

        %REJECT CHANNELS
        %[EEG, EEG.reject.delElc] = pop_rejchan(EEG,'elec',1:EEG.nbchan,'threshold',12,'norm','on','measure','kurt');

        %REJECT EPOCHS
        %EEG = pop_autorej(EEG, 'nogui','on','threshold',1000,'startprob',12);

        %REJECT CHANNELS
        %[EEG, EEG.reject.delElc] = pop_rejchan(EEG,'elec',1:EEG.nbchan,'threshold',12,'norm','on','measure','kurt');

       numcomp = numcompeig(EEG);

       EEG = pop_runica(EEG, 'extended',1,'interupt','on');%,'pca',numcomp); see https://sccn.ucsd.edu/pipermail/eeglablist/2010/003339.html
       sname = [C{1} '_' C{2} '_1st_ICA.set'];
       EEG.reject.gcompreject = ones(1,length(EEG.icachansind));
       EEG = pop_saveset(EEG,'filename',sname,'filepath',anapath); 

       clear EEG
    %end
end
