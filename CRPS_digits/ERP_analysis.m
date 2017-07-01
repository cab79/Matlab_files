clear all
close all

%grplist = [1 2 29 30]; sublist_side = {'L','R','L','R'}; sublist_grp = {'H','H','P','P'}; filepath = 'C:\Data\CRPS-DP\CRPS_Digit_Perception\';%Exp2
grplist = [35 36 37 38]; sublist_side = {'L','R','L','R'}; sublist_grp = {'H','H','P','P'}; filepath = 'C:\Data\CRPS-DP\CRPS_Digit_Perception_exp1\correcttrials\';%Exp1 left v right


run('M:\Matlab\Matlab_files\CRPS_digits\loadsubj.m');
cd(filepath)
ele_left = [];
ele_right = [];
%ele_left = [29 23 24 28 30 33 34];
%ele_right = [78 66 70 77 79 83 84];

aff_side = [2 1 1 2 1 1 1 2 1 2 1 1 1];
%aff_side = [1 1 1 1 1 1 1 1 1 1 1 1 1];

cond_nme = {'little','ring','middle','index','thumb'};
hand_nme = {'L','R'};
%sublist = [1 2];
%sublist_side = {'L','R'}; %Left, right
dp = 250;
el = 92;
latvar = 20;
no_cond = length(cond_nme); % no of conditions per data file (arm)
subjects = subjlists(grplist);
DATs = subjects;
gavgL1 = zeros(el,dp);
gavgR1 = zeros(el,dp);
countavgL1 = 0;
countavgR1 = 0;
gavgL2 = zeros(el,dp);
gavgR2 = zeros(el,dp);
countavgL2 = 0;
countavgR2 = 0;

% create TFs and grand average
sublist = {};
for s = 1:length(subjects)
    for s2 = 1:length(subjects{s,1}) 
        subj = subjects{s,1}{s2,1};
        sublist{s,s2} = subj(1:3);
        EEG = pop_loadset([subj '.set'],filepath);
        dsize = size(EEG.data);
        EEG.data = reshape(EEG.data,dsize(1),dsize(2),dsize(3));
        events = zeros(1,size(EEG.data,3));
        %j=0;
        for i = 1:size(EEG.data,3)
             if strcmp(filepath,'C:\Data\CRPS-DP\CRPS_Digit_Perception\');
                  events(1,i) = EEG.epoch(i).eventcodes{2,2};
             else
                 if length(EEG.epoch(i).eventcodes)==1 || length(EEG.epoch(i).eventcodes)==3
                    events(1,i) = EEG.epoch(i).eventcodes{1,1}{2,2};
                 elseif length(EEG.epoch(i).eventcodes)==2 || length(EEG.epoch(i).eventcodes)==4
                    events(1,i) = EEG.epoch(i).eventcodes{1,2}{2,2};
                 end
            end
        end
        cond = min(events):max(events);
        DATs{s,1}{s2,1} = cell(length(cond),1);
        %tfdata = [];
        for i = 1:no_cond
            erpdata = double(squeeze(mean(EEG.data(:,:,find(events==cond(i))),3)));
            %for e = 1:size(EEG.data,1)
            %    [tfdata(:,:,e),x,y,times,freqs] = timef(EEG.data(e,:,find(events==cond(i))),[],[-200 800],EEG.srate,1,'detret','on','winsize',60,'plotersp','off','plotitc','off','plotphase','off');   
            %end
            if strcmp(sublist_side{s},'L') && s==1
                ord = i;
                DATs{s,1}{s2,1}{ord,1} = erpdata;
                gavgL1 = gavgL1 + DATs{s,1}{s2,1}{ord,1};
                countavgL1 = countavgL1+1;
            elseif strcmp(sublist_side{s},'L') && s==3
                ord = i;
                DATs{s,1}{s2,1}{ord,1} = erpdata;
                gavgL2 = gavgL2 + DATs{s,1}{s2,1}{ord,1};
                countavgL2 = countavgL2+1;
            elseif  strcmp(sublist_side{s},'R') && s==2
                ord = (no_cond+1)-i;
                DATs{s,1}{s2,1}{ord,1} = erpdata;
                gavgR1 = gavgR1 + DATs{s,1}{s2,1}{ord,1};
                countavgR1 = countavgR1+1;
            elseif strcmp(sublist_side{s},'R') && s==4
                ord = (no_cond+1)-i;
                DATs{s,1}{s2,1}{ord,1} = erpdata;
                gavgR2 = gavgR2 + DATs{s,1}{s2,1}{ord,1};
                countavgR2 = countavgR2+1;
            end
        end
    end
end
gavgL1 = gavgL1/countavgL1;
gavgR1 = gavgR1/countavgR1;
gavgL2 = gavgL2/countavgL2;
gavgR2 = gavgR2/countavgR2;

gavg = cell(4,1);
gavg{1,1} = gavgL1;
gavg{2,1} = gavgR1;
gavg{3,1} = gavgL2;
gavg{4,1} = gavgR2;

gavg_all = (gavg{1,1}+gavg{2,1}+gavg{3,1}+gavg{4,1})/4;

save('ERP_data.mat');
no_freqs=1;
freq_idx={1};
freqs_nme = {'ERP'};
times=EEG.times;
posnegpeak = [2 1 1 2 2 2 1 1 1 1 1 1 1 1];



% normalise by common average and GFP mean over time, for each frequency
for s = 1:length(subjects)
    for f = 1:size(gavg{s,1},3)
        for e = 1:size(gavg{s,1},1)
            gavg{s,1}(e,:,f) = gavg{s,1}(e,:,f) - mean(gavg{s,1}(:,:,f),1);
        end
        gavg{s,1}(:,:,f) = gavg{s,1}(:,:,f)/mean(std(gavg{s,1}(:,:,f),1));
    end
end

% smooth gavg for peak identification
gavg_sm = gavg;
%for s = 1:length(subjects)
%    for f = 1:size(gavg{s,1},3)
%        for e = 1:size(gavg{s,1},1)
%            gavg_sm{s,1}(e,:,f) = smooth(squeeze(gavg{s,1}(e,:,f)),30,'lowess');
%        end
%    end
%end
%figure;
%imagesc(squeeze(gavg_sm{s,1}(50,:,:))'); 
%colormap jet;
%set(gca,'YDir','normal')


% identify time bins
limits=cell(no_freqs,1);
basetimesGFP = find(times<0);
stimtimesGFP = find(times>0);
timesRMSE = times(2:end);
basetimesRMSE = find(timesRMSE<0);
stimtimesRMSE = find(timesRMSE>0);

close all
for f = 1:no_freqs
    
    % GFP and CIs across subjects
    GFPg=[];
    GFP_sub=[];
    for s = 1:length(subjects)
        GFPg(:,s) = std(mean(gavg_sm{s,1}(:,:,freq_idx{f}),3),1);
        %for s2 = 1:length(subjects{s,1}) 
        %    GFP_sub(:,s,s2) = squeeze(std(squeeze(mean(subdat_all(:,:,freq_idx{f},s,s2),3)),1));
        %end
    end
    %GFP_CI=[];
    %for t = 1:length(times)
    %    GFPg(t,:) = mean(squeeze(GFP_sub(t,:,:)),2);
    %    GFP_CI(t,:) = 2.575*std(squeeze(GFP_sub(t,:,:))');
    %end
    GFPg=mean(GFPg,2);
    GFP_CI = repmat(mean(GFPg(basetimesGFP))+std(GFPg(basetimesGFP)),length(GFPg),1);
    GFP_base = repmat(mean(GFPg(basetimesGFP)),length(GFPg),1);
    figure
    yyaxis left
    plot(times,GFPg);
    hold on
    plot(times,GFP_CI,'b');
    hold on
    plot(times,GFP_base,'b--');
    
    % GFP peaks
    [peaks_GFP peaklocs_GFP] = findpeaks_Rbase(GFPg,'MinPeakHeight',unique(GFP_CI),'MinPeakProminence',unique(GFP_CI)-unique(GFP_base));
    peakloctimes_GFP = times(peaklocs_GFP);
    hold on;
    plot(peakloctimes_GFP,peaks_GFP,'v','MarkerFaceColor','b');
    % GFP troughs
    [tro_GFP trolocs_GFP] = findpeaks_Lbase(-GFPg,'MinPeakProminence',unique(GFP_CI)-unique(GFP_base));
    troloctimes_GFP = times(trolocs_GFP);
    tro_GFP = -tro_GFP;
    hold on;
    plot(troloctimes_GFP,tro_GFP,'s','MarkerFaceColor','b');
    
    % RMSE
    RMSE = [];
    RMSE_sub = [];
    for s = 1:length(subjects)
        for t = 2:size(gavg_sm{s,1},2)
            y = mean(gavg_sm{s,1}(:,t,freq_idx{f}),3);
            yhat = mean(gavg_sm{s,1}(:,t-1,freq_idx{f}),3);
            RMSE(t-1,s) = sqrt(mean((y - yhat).^2));
            %for s2 = 1:length(subjects{s,1}) 
            %    y = squeeze(mean(subdat_all(:,t,freq_idx{f},s,s2),3));
            %    yhat = squeeze(mean(subdat_all(:,t-1,freq_idx{f},s,s2),3));
            %    RMSE_sub(t-1,s,s2) = sqrt(mean((y - yhat).^2));
            %end
        end
    end
    %for t = 1:length(times)-1
    %    RMSE(t,:) = mean(squeeze(RMSE_sub(t,:,:)),2);
    %    RMSE_CI(t,:) = 2.575*std(squeeze(RMSE_sub(t,:,:))');
    %end
    RMSEav = mean(RMSE,2);
    RMSE_CI = repmat(mean(RMSEav(basetimesRMSE))+std(RMSEav(basetimesRMSE)),length(RMSEav),1);
    RMSE_base = repmat(mean(RMSEav(basetimesRMSE)),length(RMSEav),1);
    hold on
    yyaxis right
    plot(times(2:end),RMSEav,'r');
    hold on
    plot(times(2:end),RMSE_CI,'r');
    hold on
    plot(times(2:end),RMSE_base,'r--');
    
    % RMSE peaks
    [peaks_RMSE peaklocs_RMSE] = findpeaks_Rbase(RMSEav,'MinPeakHeight',unique(RMSE_CI),'MinPeakProminence',unique(RMSE_CI)-unique(RMSE_base));
    peakloctimes_RMSE = times(peaklocs_RMSE+1);
    hold on;
    plot(peakloctimes_RMSE,peaks_RMSE,'v','MarkerFaceColor','r');
    % RMSE troughs
    [tro_RMSE trolocs_RMSE] = findpeaks_Lbase(-RMSEav,'MinPeakProminence',unique(RMSE_CI)-unique(RMSE_base));
    troloctimes_RMSE = times(trolocs_RMSE+1);
    tro_RMSE = -tro_RMSE;
    hold on;
    plot(troloctimes_RMSE,tro_RMSE,'s','MarkerFaceColor','r');
    
    % find peak limits
    peakloctimes_RMSE=peakloctimes_RMSE(peakloctimes_RMSE>timesRMSE(stimtimesRMSE(1)));
    limit1=[];
    limit2=[];
    for i = 1:length(peakloctimes_RMSE)
        limit1=[limit1 peakloctimes_RMSE(i)];
        next_tros = troloctimes_RMSE(find(troloctimes_RMSE>peakloctimes_RMSE(i)));
        limit2=[limit2 next_tros(1)];
    end
    limits{f}=[limit1;limit2];
    no_peaks=size(limits{f},2);
    
    figure
    for s = 1:length(subjects)
        for n = 1:no_peaks
            subplot(length(subjects),no_peaks,(s-1)*no_peaks+n); topoplot(mean(mean(gavg{s,1}(:,find(times==limits{f}(1,n)):find(times==limits{f}(2,n)),freq_idx{f}),3),2), EEG.chanlocs);
        end
    end
end

% identify GFPs for each subject and condition
%LATgfp = NaN('double');
%GFP = NaN('double');
ugrplist = unique(sublist_grp','stable');
grplistall = reshape(repmat(ugrplist,1,size(sublist,2))',size(sublist,2)*length(ugrplist),1);
usublist = unique(sublist','stable');
usub = reshape(usublist,length(usublist)/length(ugrplist),length(ugrplist));
GFP = struct('subject',usublist,'group',grplistall);
LAT = GFP;
POW = GFP;

% lists for ANOVA LME matrix
GFP_y = [];
LAT_y = [];
POW_y = [];
G = [];
W = [];
Z = [];

avcond = {[1 5],[2 3 4]};
   
for g = 1:length(ugrplist)
    for s = 1:size(sublist,2) 
        s_name = usub{s,g};
        for f = 1:no_freqs
            f_name = freqs_nme{f};
            no_peaks = size(limits{f},2);
            for p = 1:no_peaks
                p_name = ['peak' num2str(p)];
                lim = limits{f}(:,p); % in datapoints
                for h = 1:length(hand_nme)
                    h_name = hand_nme{h};
                    for ac = 1:length(avcond)
                        c_name = ['c' num2str(ac)];
                        for i = 1:length(avcond{ac}) % for selected conditions
                            DAT(:,:,:,i) = DATs{(g-1)*2+h,1}{s,1}{avcond{ac}(i),1};
                        end
                        DAT = mean(DAT,4);

                        GFPr = std(mean(DAT(:,find(times==lim(1)):find(times==lim(2)),freq_idx{f}),3),1);
                        [~, GFPmaxLOC] = max(GFPr);
                        GFPr = sort(GFPr(:), 'descend');
                        
                        POWr = mean(DAT(:,find(times==lim(1)):find(times==lim(2)),freq_idx{f}),3);
                        if posnegpeak(p)==1
                            POWr = sort(POWr(:), 'descend');
                        elseif posnegpeak(p)==2
                            POWr = sort(POWr(:), 'ascend');
                        end

                        fieldname = [f_name '_' p_name '_' h_name '_' c_name];
                        GFP(strcmp(s_name,{GFP.subject})).(fieldname) = mean(GFPr(1:ceil(length(GFPr)*0.1)));%mean(GFPr);
                        LAT(strcmp(s_name,{LAT.subject})).(fieldname) = times(find(times==lim(1))+GFPmaxLOC-1);
                        POW(strcmp(s_name,{LAT.subject})).(fieldname) = mean(POWr(1:ceil(length(POWr)*0.1)));
                        
                        if g==1 && s==1; W = [W; f,p,h,ac];end
                    end
                end
            end
        end
    end
end

T_GFP = struct2table(GFP);
T_LAT = struct2table(LAT);
T_POW = struct2table(POW);

writetable(T_GFP,'gfp_results_erp.xlsx');
writetable(T_LAT,'lat_results_erp.xlsx');
writetable(T_POW,'pow_results_erp.xlsx');

% modify table for RM ANOVA 
T_GFP.group = categorical(T_GFP.group);
T_LAT.group = categorical(T_LAT.group);
T_POW.group = categorical(T_POW.group);
ys = {};
for ysi = 1:length(T_GFP.Properties.VariableNames(3:end))
    ys{ysi} = ['y' num2str(ysi)];
end
T_GFP.Properties.VariableNames(3:end) = ys;
T_LAT.Properties.VariableNames(3:end) = ys;
T_POW.Properties.VariableNames(3:end) = ys;


% specify and fit model
w1u = unique(W(:,1))';
for f = w1u
    w1Ind = find(W(:,1)==f);
    w2u = unique(W(w1Ind,2))';
    for p = w2u
        w2Ind = find(W(:,1)==f & W(:,2)==p);
        within = array2table(W(w2Ind,3:4));
        between = T_GFP(:,[1:2 2+w2Ind']);
        modelspec = [ys{w2Ind(1)} '-' ys{w2Ind(end)} ' ~ group'];
        GFP_rm = fitrm(between,modelspec,'WithinDesign',within);
        [ranovatbl] = ranova(GFP_rm,'WithinModel','Var1*Var2');
        if any(ranovatbl.pValue([2 5 8 11])<0.05)
            indx = find(ranovatbl.pValue([2 5 8 11])<0.05);
            sprintf('%d %d %d',f,p,indx);
        end
    end
end


