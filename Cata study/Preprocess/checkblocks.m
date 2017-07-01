%% PREPROCESSING
% 1. filtering, epoching, baseline correction, rereference
% 2. rejects chans/trials only with high frequency noise
% 3. reject ICA components related to eye movement (and 50Hz if no notch applied) but not that related to individual channel noise 
% 4. reject remaining chans/trials outside normal range
%%

clear all
filepath = 'C:\Data\Catastrophising study\Preprocessed';
cd(filepath);
files = dir('*orig.set');
load chanlocs
combine_all=0; % combining left and right stimulations, or that of different experiments, may be unwise for ICA purposes.

timebin= [-5.5 2]; % for epoching, TSOT(4)
basebin = [-5.5 -5];
stimtypes = {'c0','c1','c2','c3','c4','c5','c6','c7','c8'};
%ISIs = [1, 0.4];
filterset = [0 50];
notch_on = 1;
addpath(genpath('M:\Matlab\Matlab_files\Cata study'));
files_ana = 1:length(files);%[22,25,26,27,34,35,37,38,39];

results=cell(1,1);
count=0

for f = files_ana
    
    orig_file = files(f).name;
    [pth nme ext] = fileparts(orig_file); 
    C = strsplit(nme,'_');
    
    EEG = pop_loadset('filename',orig_file,'filepath',filepath);
    EEG.chanlocs=chanlocs;
    
  
    
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
                count=count+1;
                results{count,1} = orig_file
                results{count,2} = length(find(strcmp({EEG.event.type},'c7')))
                results{count,3} = length(find(strcmp({EEG.event.type},'c8')))
                results{count,4} = find(strcmp({EEG.event.type},'c7'));
                results{count,5} = find(strcmp({EEG.event.type},'c8'));
                %msgbox(['should be no c7 trials in block ' num2str(i)]);
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
                count=count+1;
                results{count,1} = orig_file
                results{count,2} = length(find(strcmp({EEG.event.type},'c3')))
                results{count,3} = length(find(strcmp({EEG.event.type},'c4')))
                results{count,4} = find(strcmp({EEG.event.type},'c3'));
                results{count,5} = find(strcmp({EEG.event.type},'c4'));
                %msgbox(['should be no c3 trials in block ' num2str(i)]);
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
end
save('blockcheck_results','results');