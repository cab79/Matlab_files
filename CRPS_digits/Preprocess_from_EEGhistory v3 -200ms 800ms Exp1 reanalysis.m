% 1. replicate original analysis to show this script works
% 2. need to split into left and right at the end - with a script?

clear all
if isunix
    filepath = '/scratch/cb802/Data';
    run('/scratch/cb802/Matlab_files/CRPS_digits/loadsubj.m');
else
    filepath = 'W:\Data';
    run('W:\Matlab_files\CRPS_digits\loadsubj.m');
end
ana_path1 = fullfile(filepath,'CRPS_Digit_Perception_exp1');
ana_path2 = fullfile(filepath,'CRPS_Digit_Perception');
raw_path = fullfile(filepath,'CRPS_raw/Raw');
cd(raw_path);
files = dir('*orig.set');
combine_all=0; % combining left and right stimulations, or that of different experiments, may be unwise for ICA purposes.

for f = 1:length(files)
    
    orig_file = files(f).name;
    [pth nme ext] = fileparts(orig_file); 
    C = strsplit(nme,'_');
    
    files_ana1 = dir(fullfile(ana_path1,'*.set'));
    files_ana2 = dir(fullfile(ana_path2,'*.set'));
    
    for fa = 1:length(files_ana1)
        fname = files_ana1(fa).name;
        if strfind(fname,'flip'); continue;
        elseif ~isempty(strfind(fname,C{1})) && (~isempty(strfind(fname,'Left')) || ~isempty(strfind(fname,'left')))
             EEG1l = pop_loadset(fname,ana_path1);
            EEG1l.filename = strrep(EEG1l.filename, '.left', '_left');
            EEG1l.filename = strrep(EEG1l.filename, '.Left', '_left');
            EEG1l.filename = strrep(EEG1l.filename, '.right', '_right');
            EEG1l.filename = strrep(EEG1l.filename, '.Right', '_right');
            EEG1l.filename = strrep(EEG1l.filename, '.flip', '_flip');
            EEG1l.filename = strrep(EEG1l.filename, '.Flip', '_flip');
            EEG1l.filename = strrep(EEG1l.filename, '.aff', '_aff');
            EEG1l.filename = strrep(EEG1l.filename, '.Aff', '_aff');
            EEG1l.filename = strrep(EEG1l.filename, '.Unaff', '_unaff');
            EEG1l.filename = strrep(EEG1l.filename, '.unaff', '_unaff');
            EEG1l.filename = strrep(EEG1l.filename, '_Left', '_left');
            EEG1l.filename = strrep(EEG1l.filename, '_Right', '_right');
            EEG1l.filename = strrep(EEG1l.filename, '_Flip', '_flip');
            EEG1l.filename = strrep(EEG1l.filename, '_Aff', '_aff');
            EEG1l.filename = strrep(EEG1l.filename, '_Unaff', '_unaff');
            EEG1l.filename = strrep(EEG1l.filename, '.Exp1', '_Exp1');
        elseif ~isempty(strfind(fname,C{1})) && (~isempty(strfind(fname,'Right')) || ~isempty(strfind(fname,'right')))
             EEG1r = pop_loadset(fname,ana_path1);
            EEG1r.filename = strrep(EEG1r.filename, '.left', '_left');
            EEG1r.filename = strrep(EEG1r.filename, '.Left', '_left');
            EEG1r.filename = strrep(EEG1r.filename, '.right', '_right');
            EEG1r.filename = strrep(EEG1r.filename, '.Right', '_right');
            EEG1r.filename = strrep(EEG1r.filename, '.flip', '_flip');
            EEG1r.filename = strrep(EEG1r.filename, '.Flip', '_flip');
            EEG1r.filename = strrep(EEG1r.filename, '.aff', '_aff');
            EEG1r.filename = strrep(EEG1r.filename, '.Aff', '_aff');
            EEG1r.filename = strrep(EEG1r.filename, '.Unaff', '_unaff');
            EEG1r.filename = strrep(EEG1r.filename, '.unaff', '_unaff');
            EEG1r.filename = strrep(EEG1r.filename, '_Left', '_left');
            EEG1r.filename = strrep(EEG1r.filename, '_Right', '_right');
            EEG1r.filename = strrep(EEG1r.filename, '_Flip', '_flip');
            EEG1r.filename = strrep(EEG1r.filename, '_Aff', '_aff');
            EEG1r.filename = strrep(EEG1r.filename, '_Unaff', '_unaff');
            EEG1r.filename = strrep(EEG1r.filename, '.Exp1', '_Exp1');
        end
    end
    
    for fa = 1:length(files_ana2)
        fname = files_ana2(fa).name;
        if strfind(fname,'flip'); continue;
        elseif ~isempty(strfind(fname,C{1})) && (~isempty(strfind(fname,'Left')) || ~isempty(strfind(fname,'left')))
             EEG2l = pop_loadset(fname,ana_path2);
            EEG2l.filename = strrep(EEG2l.filename, '.left', '_left');
            EEG2l.filename = strrep(EEG2l.filename, '.Left', '_left');
            EEG2l.filename = strrep(EEG2l.filename, '.right', '_right');
            EEG2l.filename = strrep(EEG2l.filename, '.Right', '_right');
            EEG2l.filename = strrep(EEG2l.filename, '.flip', '_flip');
            EEG2l.filename = strrep(EEG2l.filename, '.Flip', '_flip');
            EEG2l.filename = strrep(EEG2l.filename, '.aff', '_aff');
            EEG2l.filename = strrep(EEG2l.filename, '.Aff', '_aff');
            EEG2l.filename = strrep(EEG2l.filename, '.Unaff', '_unaff');
            EEG2l.filename = strrep(EEG2l.filename, '.unaff', '_unaff');
            EEG2l.filename = strrep(EEG2l.filename, '_Left', '_left');
            EEG2l.filename = strrep(EEG2l.filename, '_Right', '_right');
            EEG2l.filename = strrep(EEG2l.filename, '_Flip', '_flip');
            EEG2l.filename = strrep(EEG2l.filename, '_Aff', '_aff');
            EEG2l.filename = strrep(EEG2l.filename, '_Unaff', '_unaff');
            EEG2l.filename = strrep(EEG2l.filename, '.Exp1', '_Exp1');
        elseif ~isempty(strfind(fname,C{1})) && (~isempty(strfind(fname,'Right')) || ~isempty(strfind(fname,'right')))
             EEG2r = pop_loadset(fname,ana_path2);
            EEG2r.filename = strrep(EEG2r.filename, '.left', '_left');
            EEG2r.filename = strrep(EEG2r.filename, '.Left', '_left');
            EEG2r.filename = strrep(EEG2r.filename, '.right', '_right');
            EEG2r.filename = strrep(EEG2r.filename, '.Right', '_right');
            EEG2r.filename = strrep(EEG2r.filename, '.flip', '_flip');
            EEG2r.filename = strrep(EEG2r.filename, '.Flip', '_flip');
            EEG2r.filename = strrep(EEG2r.filename, '.aff', '_aff');
            EEG2r.filename = strrep(EEG2r.filename, '.Aff', '_aff');
            EEG2r.filename = strrep(EEG2r.filename, '.Unaff', '_unaff');
            EEG2r.filename = strrep(EEG2r.filename, '.unaff', '_unaff');
            EEG2r.filename = strrep(EEG2r.filename, '_Left', '_left');
            EEG2r.filename = strrep(EEG2r.filename, '_Right', '_right');
            EEG2r.filename = strrep(EEG2r.filename, '_Flip', '_flip');
            EEG2r.filename = strrep(EEG2r.filename, '_Aff', '_aff');
            EEG2r.filename = strrep(EEG2r.filename, '_Unaff', '_unaff');
            EEG2r.filename = strrep(EEG2r.filename, '.Exp1', '_Exp1');
        end
    end
    
    hist=cell(4,1);
    hist{1,1} = EEG1l.history;
    hist{2,1} = EEG1r.history;
    hist{3,1} = EEG2l.history;
    hist{4,1} = EEG2r.history;
    
    rejcall=[];
    timeall=cell(4,1);
    for h = 1:length(hist)
        his = regexp(hist{h},'([^;]*)','tokens');
        his = [his{:}];
        for i=1:length(his)
            if ~isempty(strfind(his{i},'pop_repchan'))
                rejc = regexp(his{i},'\[(.*?)\]','tokens');
                if isempty(rejc)
                    rejc = strsplit(his{i},',');
                    rejc = rejc(3);
                else
                    rejc = rejc{:};
                end
                rejcall = [rejcall str2num(rejc{:})];
            elseif ~isempty(strfind(his{i},'time')) && ~isempty(strfind(his{i},'pop_select'))
                timeall{h} = regexp(his{i},'\[(.*?)\]','tokens');
            end
        end
    end
    
    
    EEG = pop_loadset('filename',orig_file,'filepath',raw_path);
    EEG = pop_eegfilt( EEG, 0.5, 0, [], 0);
    EEG = pop_eegfilt( EEG, 0, 30, [], 0);
    EEG = pop_eegfilt(EEG,49,51,[],1);
    
    
    if ~isempty(rejcall); 
        try
            EEG = eeg_interp(EEG,unique(rejcall),'spherical');
        end
    end
    
    %try
        [EEG, EEG.reject.delElc] = pop_rejchan(EEG,'elec',1:92,'threshold',5,'norm','on','measure','kurt');
    %    EEG = eeg_interp(EEG,EEG.reject.delElc,'spherical');
    %catch
    %    try
    %        [EEG, EEG.reject.delElc] = pop_rejchan(EEG,'elec',1:92,'threshold',6,'norm','on','measure','kurt');
    %        EEG = eeg_interp(EEG,EEG.reject.delElc,'spherical');
    %        try
    %            [EEG, EEG.reject.delElc] = pop_rejchan(EEG,'elec',1:92,'threshold',7,'norm','on','measure','kurt');
    %            EEG = eeg_interp(EEG,EEG.reject.delElc,'spherical');
    %        end
    %    end
    %end
    
    time_exp1 = timeall{1};
    %time_exp2 = timeall{3};
    time_exp1 = time_exp1{:};
    %time_exp2 = time_exp2{:};
    EEGn1 = pop_select( EEG,'time',str2num(time_exp1{:}));
    %EEGn2 = pop_select( EEG,'time',str2num(time_exp2{:}));
    
    %EEG=EEGn1;
    %part1accu;
    %EEGn1=EEG;
    
    EEGn1 = pop_epoch( EEGn1, {'STIM'}, [-0.2 0.8], 'newname', [C{1} '_epochs1'], 'epochinfo', 'yes');
    %EEGn2 = pop_epoch( EEGn2, {'STIM'}, [-0.2 0.8], 'newname', [C{1} '_epochs2'], 'epochinfo', 'yes');
    EEGn1 = pop_rmbase( EEGn1, [-200    0]);
    %EEGn2 = pop_rmbase( EEGn2, [-200    0]);
    
    EEGn1 = pop_reref( EEGn1, []);
    EEGn1 = pop_autorej(EEGn1, 'nogui','on','threshold',100,'startprob',5);
    
    %exp1
   %allE=[];
   %for ep = 1:length(EEGn1.epoch)
   %     stimevidx = find(strcmp('STIM',EEGn1.epoch(ep).eventtype));
   %     if ~isempty(stimevidx)
   %         allE(ep) = EEGn1.epoch(ep).eventinit_index{stimevidx};
   %     end
   %end
   %keptEl=[];
   %for ep = 1:length(EEG1l.epoch)
   %     stimevidx = find(strcmp('STIM',EEG1l.epoch(ep).eventtype));
   %     if ~isempty(stimevidx)
   %         keptEl(ep) = EEG1l.epoch(ep).eventinit_index{stimevidx};
   %     end
   %end
   %keptEr=[];
   %for ep = 1:length(EEG1r.epoch)
   %     stimevidx = find(strcmp('STIM',EEG1r.epoch(ep).eventtype));
   %     if ~isempty(stimevidx)
   %         keptEr(ep) = EEG1r.epoch(ep).eventinit_index{stimevidx};
   %     end
   %end
   %keptEl(keptEl==0)=[];
   %keptEr(keptEr==0)=[];
   %if combine_all == 1;
   %    keptE = sort([keptEl keptEr]);
   %    [inE iall ikept] = intersect(allE,keptE);
   %    EEGn1 = pop_select(EEGn1,'trial',iall);
   %    n1_trials = EEGn1.trials;
   %else
   %    [inE iall ikept] = intersect(allE,keptEl);
   %    EEGn1l = pop_select(EEGn1,'trial',iall);
   %    n1l_trials = EEGn1l.trials;
   %    [inE iall ikept] = intersect(allE,keptEr);
   %    EEGn1r = pop_select(EEGn1,'trial',iall);
   %    n1r_trials = EEGn1r.trials;
   %end
   
   %exp2
   %allE = [EEGn2.epoch.eventinit_index];
   %keptEl = [EEG2l.epoch.eventinit_index];
   %keptEr = [EEG2r.epoch.eventinit_index];
   %if combine_all == 1;
   %    keptE = sort([keptEl keptEr]);
   %    [inE iall ikept] = intersect(allE,keptE);
   %    EEGn2 = pop_select(EEGn2,'trial',iall);
   %    n2_trials = EEGn2.trials;
   %else
   %    [inE iall ikept] = intersect(allE,keptEl);
   %    EEGn2l = pop_select(EEGn2,'trial',iall);
   %    n2l_trials = EEGn2l.trials;
   %    [inE iall ikept] = intersect(allE,keptEr);
   %    EEGn2r = pop_select(EEGn2,'trial',iall);
   %    n2r_trials = EEGn2r.trials;
   %end
   
   if combine_all == 1;
       EEG = pop_mergeset(EEGn1, EEGn2, 0);
       %EEG = pop_reref( EEG, []);
       %correctRank = EEG.nbchan-1;
       %EEG = pop_runica(EEG, 'extended',1,'interupt','on','pca',20);
       sname = [C{1} '_ICA_30Hz.set'];
       EEG = pop_saveset(EEG,'filename',sname,'filepath',raw_path); 
       save([C{1} '_trials_30Hz.mat'],'trials');
   else 
       %EEG = pop_reref( EEGn1l, []);
       %correctRank = EEG.nbchan-1;
       %EEG = pop_runica(EEG, 'extended',1,'interupt','on','pca',20);
       EEG=EEGn1;
       sname = [C{1} '_Exp1.set'];
       EEG = pop_saveset(EEG,'filename',sname,'filepath',raw_path); 
       %trials = n1l_trials;
       %save([C{1} '_Exp1_left_trials_30Hz.mat'],'trials');
       EEGc=EEG;
         selectfnum =1:5;
        EEG=EEGc;
        part1analysis
       correctRank = EEG.nbchan-1;
       EEG = pop_runica(EEG, 'extended',1,'interupt','on','pca',20);
       sname = [C{1} '_Exp1_left.set'];
       EEG = pop_saveset(EEG,'filename',sname,'filepath',raw_path); 
       %trials = n1l_trials;
       %save([C{1} '_Exp1_left_trials_30Hz.mat'],'trials');

       selectfnum =6:10;
       EEG=EEGc;
        part1analysis
       correctRank = EEG.nbchan-1;
       EEG = pop_runica(EEG, 'extended',1,'interupt','on','pca',20);
       sname = [C{1} '_Exp1_right.set'];
       EEG = pop_saveset(EEG,'filename',sname,'filepath',raw_path); 
       
       %EEG = pop_reref( EEGn1r, []);
       %correctRank = EEG.nbchan-1;
       %EEG = pop_runica(EEG, 'extended',1,'interupt','on','pca',20);
       %sname = [C{1} '_Exp1_right.set'];
       %EEG = pop_saveset(EEG,'filename',sname,'filepath',raw_path); 
       %trials = n1r_trials;
       %save([C{1} '_Exp1_right_trials_30Hz.mat'],'trials');
       
   %    EEG = pop_reref( EEGn2l, []);
   %    correctRank = EEG.nbchan-1;
   %    EEG = pop_runica(EEG, 'extended',1,'interupt','on','pca',20);
   %    sname = [C{1} '_Exp2_left_ICA_30Hz.set'];
   %    EEG = pop_saveset(EEG,'filename',sname,'filepath',raw_path); 
   %    trials = n2l_trials;
   %    save([C{1} '_Exp2_left_trials_30Hz.mat'],'trials');
       
   %    EEG = pop_reref( EEGn2r, []);
   %    correctRank = EEG.nbchan-1;
   %    EEG = pop_runica(EEG, 'extended',1,'interupt','on','pca',20);
   %    sname = [C{1} '_Exp2_right_ICA_30Hz.set'];
   %    EEG = pop_saveset(EEG,'filename',sname,'filepath',raw_path); 
   %    trials = n2r_trials;
   %    save([C{1} '_Exp2_right_trials_30Hz.mat'],'trials');
   end
   
   clear EEG EEG2r EEG2l EEG1r EEG1l EEGn1l EEGn1r EEGn2l EEGn2r
       
   
   
   %EEG = pop_reref( EEG, []);
   %split into exp1 and exp2
   
   %EEG12=EEG;
   %selectfnum =1:5;
   %part2analysis;
   %EEG2l=EEG;
   
   %EEG=EEG12;
   %selectfnum =6:10;
   %part2analysis;
   %EEG2r=EEG;
   
   
   
   
   
   
   
   
   
   
   %icaweights = EEG2l.icaweights;
   %icasphere = EEG2l.icasphere;
   %data2l = reshape(EEG2l.data,EEG2l.nbchan,EEG2l.trials*EEG2l.pnts);
   %icaactL = (icaweights)*data2l(EEG2l.icachansind,:);
   
   %icaweights = EEG2r.icaweights;
   %icasphere = EEG2r.icasphere;
   %data2r = reshape(EEG2r.data,EEG2r.nbchan,EEG2r.trials*EEG2r.pnts);
   %icaactR = (icaweights)*data2r(EEG2r.icachansind,:);
   
   %ICAw = EEG2l.icaweights;
   %ICAs = EEG2l.icasphere;
   %dataL = reshape(EEGl.data,EEGl.nbchan,EEGl.trials*EEGl.pnts);
   %ICAactL = (ICAw)*dataL(EEG2l.icachansind,:);
   
   %ICAw = EEG2r.icaweights;
   %ICAs = EEG2r.icasphere;
   %dataR = reshape(EEGr.data,EEGr.nbchan,EEGr.trials*EEGr.pnts);
   %ICAactR = (ICAw)*dataR(EEG2r.icachansind,:);
   
   %cormat = [];
   %thresh=[0.7];
   %datcorsL=[];
   %datcorsR=[];
   %for t=1:length(thresh)
   %    for i = 1:size(icaactL,1)
   %        cormat(i) = corr(ICAactL(i,:)',icaactL(i,:)','type','Spearman');
   %    end
   %    retain=find(cormat>thresh(t));
   %    length(retain)%%

   %    rej = 1:size(icaactL,1);
   %    rej(retain)=[];
   %    EEG2l_t = pop_subcomp( EEG2l, rej, 0);
   %    EEG2r_t = pop_subcomp( EEG2r, rej, 0);
       
   %    dataL = reshape(EEGl.data,EEGl.nbchan,EEGl.trials*EEGl.pnts);
   %    dataR = reshape(EEGr.data,EEGr.nbchan,EEGr.trials*EEGr.pnts);
   %    dataLn = reshape(EEG2l_t.data,EEG2l_t.nbchan,EEG2l_t.trials*EEG2l_t.pnts);
   %    dataRn = reshape(EEG2r_t.data,EEG2r_t.nbchan,EEG2r_t.trials*EEG2r_t.pnts);
       
   %    datcorL=[];
   %    datcorR=[];
   %    for i=1:EEGl.nbchan
   %        datcorL(i) = corr(dataL(i,:)',dataLn(i,:)');
   %        datcorR(i) = corr(dataR(i,:)',dataRn(i,:)');
   %    end
   %    datcorsL(t) = mean(datcorL)
   %    datcorsR(t) = mean(datcorR)
       
   %end
   
   
end
matlabmail
