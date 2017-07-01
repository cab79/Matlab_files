%% PREPROCESSING
% 1. filtering, epoching, baseline correction, rereference
% 2. rejects chans/trials only with high frequency noise
% 3. reject ICA components related to eye movement (and 50Hz if no notch applied) but not that related to individual channel noise 
% 4. reject remaining chans/trials outside normal range
%%

clear all
filepath = 'C:\Data\Catastrophising study\Orig';
anapath = 'C:\Data\Catastrophising study\Preprocessed';
cd(filepath);
files = dir('*orig.set');
load('C:\Data\Catastrophising study\Orig\chanlocs.mat')
combine_all=0; % combining left and right stimulations, or that of different experiments, may be unwise for ICA purposes.

timebin= [-5.5 2]; % for epoching, TSOT(4)
basebin = [-5.5 -5];
stimtypes = {'c0','c1','c2','c3','c4','c5','c6','c7','c8'};
%ISIs = [1, 0.4];
filterset = [0 50];
notch_on = 1;
detrend_on = 1;
addpath(genpath('M:\Matlab\Matlab_files\Cata study'));
files_ana = [1]%length(files);%[22,25,26,27,34,35,37,38,39];


for f = files_ana
    
    orig_file = files(f).name;
    [pth nme ext] = fileparts(orig_file); 
    C = strsplit(nme,'_');
    
    EEG = pop_loadset('filename',orig_file,'filepath',filepath);
    EEG.chanlocs=chanlocs;
    
    % FILTER
    if notch_on==1 
        if filterset(2)>45; EEG = pop_eegfiltnew(EEG,45,55,[],1); end
        if filterset(2)>95; EEG = pop_eegfiltnew(EEG,95,105,[],1); end
    elseif notch_on==2
        EEG.data = rm50Hz(EEG.data,[2 1 3],EEG.srate,20);
        %pause
        close all
    end;
    if filterset(1)>0; EEG = pop_eegfiltnew( EEG, filterset(1), 0, [], 0);end
    if filterset(2)>0; EEG = pop_eegfiltnew( EEG, 0, filterset(2), [], 0);end
    
    % replace event types (see "Conditions.xls")
    for i = 1:length(EEG.event)
        if EEG.event(i).type=='1'
            EEG.event(i).type='c0';
        elseif EEG.event(i).type=='3'
            EEG.event(i).type='c3';
        elseif EEG.event(i).type=='5'
            EEG.event(i).type='c4';
        elseif EEG.event(i).type=='6'
            EEG.event(i).type='c7';
        elseif EEG.event(i).type=='7'
            EEG.event(i).type='c8';
        end
    end
    
    %--separate into blocks to identify further conditions--%
    
    %option 1: boundary
    %boun = strcmp({EEG.event.type},'boundary');
    %for i = 1:length(boun)-1
    %    if boun(i)==1 && boun(i+1)==1
    %        boun(i)=0;
    %    end
    %end
    
    %option 2: large gaps
    lat = [EEG.event.latency]; 
    latdiff = lat(2:end)-lat(1:end-1);
    boun = latdiff>10000;
    
    fboun = find(boun);
    fboun = [fboun length(EEG.event)];
    
    EEGall=EEG;
    newEEGall=[];
    for i = 1:length(fboun)-1
        EEG = pop_select(EEGall,'point',[EEGall.event(fboun(i)).latency-timebin(1)*EEGall.srate EEGall.event(fboun(i+1)).latency+timebin(2)*EEGall.srate]);
        if (length(find(strcmp({EEG.event.type},'c3')))+length(find(strcmp({EEG.event.type},'c4')))) + (length(find(strcmp({EEG.event.type},'c7')))+length(find(strcmp({EEG.event.type},'c8'))))==0
            continue
        end
        
        if (length(find(strcmp({EEG.event.type},'c3')))+length(find(strcmp({EEG.event.type},'c4')))) > (length(find(strcmp({EEG.event.type},'c7')))+length(find(strcmp({EEG.event.type},'c8')))) % without-task block
            if length(find(strcmp({EEG.event.type},'c7')))+length(find(strcmp({EEG.event.type},'c8')))>0
                msgbox(['should be no c7 trials in block ' num2str(i)]);
            end
            for ii = 1:length(EEG.event)
                if EEG.event(ii).type=='2'
                    EEG.event(ii).type='c1';
                elseif EEG.event(ii).type=='4'
                    EEG.event(ii).type='c2';
                end
            end
        elseif (length(find(strcmp({EEG.event.type},'c3')))+length(find(strcmp({EEG.event.type},'c4')))) < (length(find(strcmp({EEG.event.type},'c7')))+length(find(strcmp({EEG.event.type},'c8')))) % with-task block
            if length(find(strcmp({EEG.event.type},'c3')))+length(find(strcmp({EEG.event.type},'c4')))>0
                msgbox(['should be no c3 trials in block ' num2str(i)]);
            end
            for ii = 1:length(EEG.event)
                if EEG.event(ii).type=='2'
                    EEG.event(ii).type='c5';
                elseif EEG.event(ii).type=='4'
                    EEG.event(ii).type='c6';
                end
            end
        end
        if i==1 || isempty(newEEGall)
            newEEGall = EEG;
        else
            newEEGall = pop_mergeset(newEEGall, EEG);
        end
    end
    EEG = newEEGall;
        
    % EPOCH
    EEG = pop_epoch( EEG, stimtypes, timebin, 'newname', [C{1} '_' C{2} '_epochs'], 'epochinfo', 'yes');
        
    % RE-REFERENCE
    EEG = pop_reref( EEG, []);
    
    % LINEAR DETREND
    if detrend_on
        tempdata = [];
        %trend = [];
        for i = 1:EEG.trials
            tempdata(:,:,i) = detrend(EEG.data(:,:,i)')'; 
            %trend(:,:,i) = EEG.data(:,:,i)-tempdata(:,:,i);
            EEG.data(:,:,i) = tempdata(:,:,i);
        end
        clear tempdata
    end
    
    % REMOVE BASELINE
    EEG = pop_rmbase( EEG, [basebin(1)*1000 basebin(2)*1000]);
    
    sname = [C{1} '_' C{2} '_epoched.set'];
    EEG = pop_saveset(EEG,'filename',sname,'filepath',anapath); 
   
    %orig_epochs = [EEG.epoch.eventurevent];
    %save(fullfile(anapath,[C{1} '_trend.mat']),'trend','orig_epochs');
end


for f = files_ana
    orig_file = files(f).name;
    [pth nme ext] = fileparts(orig_file); 
    C = strsplit(nme,'_');
    lname = [C{1} '_' C{2} '_epoched.set'];
    EEG = pop_loadset('filename',lname,'filepath',anapath);
    
    EEG = FTrejman(EEG,[]); % high freq to identify noise not related to eye movement
    
    sname = [C{1} '_' C{2} '_manrej.set'];
    EEG = pop_saveset(EEG,'filename',sname,'filepath',anapath); 
end


for f = files_ana
    
    orig_file = files(f).name;
    [pth nme ext] = fileparts(orig_file); 
    C = strsplit(nme,'_');
    lname = [C{1} '_' C{2} '_manrej.set'];
    EEG = pop_loadset('filename',lname,'filepath',anapath);
    
    %REJECT CHANNELS
   % [EEG, EEG.reject.delElc] = pop_rejchan(EEG,'elec',1:EEG.nbchan,'threshold',12,'norm','on','measure','kurt');
    
    %REJECT EPOCHS
   % EEG = pop_autorej(EEG, 'nogui','on','threshold',1000,'startprob',12);
    
    %REJECT CHANNELS
   % [EEG, EEG.reject.delElc] = pop_rejchan(EEG,'elec',1:EEG.nbchan,'threshold',12,'norm','on','measure','kurt');
   
   numcomp = numcompeig(EEG);
   EEG = pop_runica(EEG, 'extended',1,'interupt','off','pca',numcomp); %see https://sccn.ucsd.edu/pipermail/eeglablist/2010/003339.html
   sname = [C{1} '_' C{2} '_ICA.set'];
   EEG = pop_saveset(EEG,'filename',sname,'filepath',anapath); 

   clear EEG
end


for f = files_ana

    orig_file = files(f).name;
    [pth nme ext] = fileparts(orig_file); 
    C = strsplit(nme,'_');
    
    lname = [C{1} '_' C{2} '_ICA.set'];
    EEG = pop_loadset('filename',lname,'filepath',anapath)
        
       
   EEG = autorejcomp(EEG,1);
   [EEG temp detrendcomps] = autorejcomp(EEG,6,[0 0.8],1,[],[-0.5 0],3,0, 0.30);

   EEG = pop_saveset(EEG,'filename',lname,'filepath',anapath); 
   save(fullfile(anapath,[C{1} '_detrendcomps']),'detrendcomps');

   clear EEG
end
eeglab
