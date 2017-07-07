%% PREPROCESSING
% 1. filtering, epoching, baseline correction, rereference
% 2. rejects chans/trials only with high frequency noise
% 3. reject ICA components related to eye movement (and 50Hz if no notch applied) but not that related to individual channel noise 
% 4. reject remaining chans/trials outside normal range
%%

clear all
origpath = 'C:\Data\CORE\Orig\';
anapath = 'C:\Data\CORE\Preprocessed_100Hz';
design_path = 'C:\Data\CORE\Raw';
cd(origpath);
files = dir('*orig.set');
cd(anapath);
combine_all=0; % combining left and right stimulations, or that of different experiments, may be unwise for ICA purposes.

timebins = [-0.2 0.9; % for epoching, TSOT(2)
            -0.2 0.3]; % for epoching, TSOT(4)
ISIs = [1, 0.4];
filterset = [0.5 100];
%Sr = 500;
notch_on = 1;
addpath(genpath('C:\Matlab_files\CORE\Supporting_functions'));

markers = {'STIM' 'DIN2'};
use_marker = 2;
STIM_delay = 0;
leftmarkers = [1:4 9:12 17:20]; 
rightmarkers = [5:8 13:16 21:24]; 
targets = [1 2 5 6 9 10 13 14 17 18 21 22];
nontargets = [3 4 7 8 11 12 15 16 19 20 23 24];

files_ana = [68 69]%1:length(files);
trials_ana = 1; fname_ext = '';

for f = files_ana
    
    orig_file = files(f).name;
    [pth nme ext] = fileparts(orig_file); 
    C = strsplit(nme,'_');
    
    EEG = pop_loadset('filename',orig_file,'filepath',origpath);
    
    if strcmp(C{2},'2')
        timebin = timebins(1,:);
        isi = ISIs(1);    
        %EEG=pop_resample(EEG,Sr);
    elseif strcmp(C{2},'4')
        timebin = timebins(2,:);
        isi = ISIs(2);
    end
    
    if strcmp(C{1},'CORE000') || strcmp(C{1},'CORE001') || strcmp(C{1},'CORE002')
        for e = 1:length(EEG.event)
            if strcmp(EEG.event(e).type,'DIN2')
                EEG.event(e).latency = EEG.event(e).latency-112;
            end
        end
    end
    
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
    %if strcmp(C{1},'CORE020')
    %    desfile = dir(fullfile(design_path,C{1}, ['*part' C{2} '*']));
    %    load(fullfile(design_path,C{1},desfile.name));
    %    des = dt.design;
    %end
    if strcmp(markers{use_marker},'DIN2')
        rm_trials = [];
        for ep = 1:length(EEG.epoch)
            stimevidx = find(strcmp('STIM',EEG.epoch(ep).eventtype));
            if ep<length(EEG.epoch); stimevidx1 = find(strcmp('STIM',EEG.epoch(ep+1).eventtype));end;
            if isempty(stimevidx)
                if length(stimevidx1)>1
                    event = EEG.epoch(ep+1).eventurevent{1, stimevidx1(1)};
                    event_idx = find([EEG.event.urevent]==event);
                    %EEG = pop_editeventvals(EEG, 'changefield', {event_idx 'latency' eeg_point2lat(EEG.event(event_idx-1).latency,[],EEG.srate, [EEG.xmin EEG.xmax])});
                    EEG.event(event_idx).latency = EEG.event(event_idx-1).latency;
                else
                    rm_trials = [rm_trials ep];
                %
                %    try 
                        
                %    catch
                %        error(['not enough stims in epoch ' num2str(ep+1)]);
                %    end
                end
            end
        end
        EEG = pop_select(EEG,'notrial',rm_trials);
        % EPOCH AGAIN AFTER STIM CORRECTION
        EEG = pop_epoch( EEG, {'DIN2'}, timebin, 'newname', [C{1} '_' C{2} '_epochs'], 'epochinfo', 'yes');
    end


    EEG = rmTrailsISIvar(EEG,0.003,timebin,isi,markers{use_marker});
        
    % RE-REFERENCE
    EEG = pop_reref( EEG, []);
    
    % LINEAR DETREND
    for i = 1:EEG.trials, EEG.data(:,:,i) = detrend(EEG.data(:,:,i)')'; end;
    
    % REMOVE BASELINE
    EEG = pop_rmbase( EEG, [timebin(1)*100    0]);
    
    sname = [C{1} '_' C{2} '_epoched.set'];
    EEG = pop_saveset(EEG,'filename',sname,'filepath',anapath); 
    
end

%remove = struct;
%remove.scales = 1:3;
%remove.times = [-0.03 0.03];
%keep = struct;
%keep.scales = 1:3;
%keep.times = [-0.03 0.03];
%EEG = EP_den_EEG(EEG,5,[],'create_den_coeff',remove,[],34);

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
    
    if strcmp(C{2},'2')
        timebin = timebins(1,:);
    elseif strcmp(C{2},'4')
        timebin = timebins(2,:);
    end
    
    [conds, tnums, fnums, bnums] = get_markers(EEG);
    
    EEGall = EEG;
    for hand = 1:2
        
        if hand == 1
            Hmarkers = leftmarkers;
            handname = 'left';
            elec = 70;
        else
            Hmarkers = rightmarkers;
            handname = 'right';
            elec = 34;
        end
        
        % create an index of conds
        selectepochs=[];
        for i = 1:length(Hmarkers)
            se = find(conds==Hmarkers(i));
            selectepochs = [selectepochs se(1:ceil(length(se)*trials_ana))];
        end
        selectepochs = sort(selectepochs);
    
        EEG = pop_select(EEGall,'trial',selectepochs);
        
        %for t = 1:24
        %    no_trials = ceil(length(EEG.epoch)*trials_ana);
        %    selecttrials = 1:no_trials;
        %    EEG = pop_select(EEG,'trial',selecttrials);
        %end

        %REJECT CHANNELS
        %[EEG, EEG.reject.delElc] = pop_rejchan(EEG,'elec',1:EEG.nbchan,'threshold',12,'norm','on','measure','kurt');

        %REJECT EPOCHS
        %EEG = pop_autorej(EEG, 'nogui','on','threshold',1000,'startprob',12);

        %REJECT CHANNELS
        %[EEG, EEG.reject.delElc] = pop_rejchan(EEG,'elec',1:EEG.nbchan,'threshold',12,'norm','on','measure','kurt');
        
        
        if strcmp(C{2},'2')
            close all
            markind = find(strcmp({EEG.event.type},markers(use_marker)));
            [Hconds, Htnums, Hfnums, Hbnums] = get_markers(EEG);
            event_targs = [];
            for t = 1:length(targets)
                event_targs = [event_targs find(Hconds==targets(t))];
            end
            event_nontargs = [];
            for t = 1:length(nontargets)
                event_nontargs = [event_nontargs find(Hconds==nontargets(t))];
            end
            event_targs = sort(event_targs);
            event_nontargs = sort(event_nontargs);
            %EEG = acstp_on_epochs(EEG,elec,markers(use_marker),[],[0.03,0.2],0,[],0); % no split of targs and non-targs
            EEG = acstp_on_epochs(EEG,elec,{markind(event_targs)},{markind(event_nontargs)},[],[0.03,0.5],0,[],1); % split of targs and non-targs
            sname = [C{1} '_' C{2} '_' handname fname_ext '_ACSTP.set'];
            EEG = pop_saveset(EEG,'filename',sname,'filepath',anapath);
        elseif strcmp(C{2},'4')
           numcomp = numcompeig(EEG);
           EEG = pop_runica(EEG, 'extended',1,'interupt','on');%,'pca',numcomp); see https://sccn.ucsd.edu/pipermail/eeglablist/2010/003339.html
           sname = [C{1} '_' C{2} '_' handname fname_ext '_1st_ICA.set'];
           EEG = pop_saveset(EEG,'filename',sname,'filepath',anapath); 
        end
        
       clear EEG
    end
end


for f = files_ana

    orig_file = files(f).name;
    [pth nme ext] = fileparts(orig_file); 
    C = strsplit(nme,'_');
    
    if strcmp(C{2},'2')
        timebin = timebins(1,:);
        continue
    elseif strcmp(C{2},'4')
        timebin = timebins(2,:);
    end

    for hand = 1:2
        
        if hand == 1
            handname = 'left';
        else
            handname = 'right';
        end

        lname = [C{1} '_' C{2} '_' handname fname_ext '_1st_ICA.set'];
        EEG = pop_loadset('filename',lname,'filepath',anapath);
        
       
       EEG = autorejcomp(EEG,1);
       [EEG temp] = autorejcomp(EEG,6,[0.03 0.23],1,45,[timebin(1) 0],3,0, 0.30);
       
       EEG = pop_saveset(EEG,'filename',lname,'filepath',anapath); 

       clear EEG
    end
end
eeglab
