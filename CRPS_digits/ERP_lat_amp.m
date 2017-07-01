clear all
close all
%grplist = [1 3 10 12]; sublist_side = {'L','R','L','R'}; %flipped
%grplist = [33 34 31 32]; sublist_side = {'L','R','L','R'}; %Left, right, OR mostly left, mostly right %unflipped
%grplist = [1 2 29 30]; sublist_side = {'L','R','L','R'}; 
grplist = [35 36 37 38]; sublist_side = {'L','R','L','R'}; %Exp1 left v right
%grplist = [43 44];% 45 46]; 
%grplist = [47:50]; %Exp2 left v right - NEW ANALYSIS
%sublist_side = {'L','R','L','R'}; 
%grplist = [37 38]; sublist_side = {'L','R'}
%grplist = [39 40 41 42]; sublist_side = {'L','R','L','R'}; %Affected vs unaffected exp1
%filepath = 'W:\Data\CRPS_Digit_Perception\';
filepath = 'W:\Data\CRPS_Digit_Perception_exp1\alltrials';
%ele_left = [29 23 24 28 30 33 34];
%ele_right = [78 66 70 77 79 83 84];
ele_left = [];
ele_right = [];
%aff_side = [2 1 1 2 1 1 1 2 1 2 1 1 1];
aff_side = [1 1 1 1 1 1 1 1 1 1 1 1 1];
flip=0;

peaks_nme = {'1','2','3','4','5','6','7'};
%peaks_wind = {20,48;
%              48,88;
%              88,132;
%              168,228;
%              280,340;
%              552,748;};
peaks_wind = {20,48;
              48,88;
              88,120;
              152,200;
              252,352; % for exp2: 200,352; for exp1: 252,352;
              552,748;}; 
cond_nme = {'little','ring','middle','index','thumb'};
elegrp_nme = {'Lpos','Rpos','Lneg','Rneg'};
hand_nme = unique(sublist_side, 'stable');

%sublist = [1 2];
%sublist_side = {'L','R'}; %Left, right
dp = 250;
el = 92;
no_ele = 5; % number of electrodes to use to identify latencies
no_cond = length(cond_nme); % no of conditions per data file (arm)
no_peaks = length(peaks_nme);
no_elegrp = length(elegrp_nme);
end_search = 360; % time (ms) to use as end limit for final peak 


run('W:\Matlab_files\CRPS_digits\loadsubj');
subjects = subjlists(grplist);
ERPs = subjects;
gavgL1 = zeros(el,dp);
gavgR1 = zeros(el,dp);
countavgL1 = 0;
countavgR1 = 0;
gavgL2 = zeros(el,dp);
gavgR2 = zeros(el,dp);
countavgL2 = 0;
countavgR2 = 0;

% create ERPs and grand average
sublist = {};
for s = 1:length(subjects)
    for s2 = 1:length(subjects{s,1}) 
        subj = subjects{s,1}{s2,1};
        sublist{s,s2} = subj(1:3);
        
        tmp_nme = subj;
          tmp_nme = strrep(tmp_nme, '.left', '_left');
            tmp_nme = strrep(tmp_nme, '.Left', '_left');
            tmp_nme = strrep(tmp_nme, '.right', '_right');
            tmp_nme = strrep(tmp_nme, '.Right', '_right');
            tmp_nme = strrep(tmp_nme, '.flip', '_flip');
            tmp_nme = strrep(tmp_nme, '.Flip', '_flip');
            tmp_nme = strrep(tmp_nme, '.aff', '_aff');
            tmp_nme = strrep(tmp_nme, '.Aff', '_aff');
            tmp_nme = strrep(tmp_nme, '.Unaff', '_unaff');
            tmp_nme = strrep(tmp_nme, '.unaff', '_unaff');
            tmp_nme = strrep(tmp_nme, '_Left', '_left');
            tmp_nme = strrep(tmp_nme, '_Right', '_right');
            tmp_nme = strrep(tmp_nme, '_Flip', '_flip');
            tmp_nme = strrep(tmp_nme, '_Aff', '_aff');
            tmp_nme = strrep(tmp_nme, '_Unaff', '_unaff');
            tmp_nme = strrep(tmp_nme, '.Exp1', '_Exp1');
            tmp_nme = strrep(tmp_nme, '.Exp2', '_Exp2');
        subj=tmp_nme;
        EEG = pop_loadset([subj '.set'],filepath);
        if flip==1 && (~isempty(strfind(subj, 'right')) || ~isempty(strfind(subj, 'Right'))); EEG = flipchan(EEG); end;
        dsize = size(EEG.data);
        %EEG = pop_eegfilt(EEG,0,25,0,0);
        EEG.data = reshape(EEG.data,dsize(1),dsize(2),dsize(3));
        %EEG.data = iirfilt(EEG.data,250,0,25);
        events = zeros(1,size(EEG.data,3));
        %j=0;
        for i = 1:size(EEG.data,3)
            %if strcmp(EEG.event(i).type,'STIM')
                %j = j+1;
             if strcmp(filepath,'W:\Data\CRPS_Digit_Perception\') || ~isempty(strfind(subj,'Exp2'));
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
        ERPs{s,1}{s2,1} = cell(length(cond),1);
        chanlocs = EEG.chanlocs;
        for i = 1:no_cond
            if strcmp(sublist_side{s},'L') && s==1
                ord = i;
                ERPs{s,1}{s2,1}{ord,1} = double(squeeze(mean(EEG.data(:,:,find(events==cond(i))),3))); % list, subject, condition
                gavgL1 = gavgL1 + ERPs{s,1}{s2,1}{ord,1};
                countavgL1 = countavgL1+1;
            elseif strcmp(sublist_side{s},'L') && s==3
                ord = i;
                ERPs{s,1}{s2,1}{ord,1} = double(squeeze(mean(EEG.data(:,:,find(events==cond(i))),3))); % list, subject, condition
                gavgL2 = gavgL2 + ERPs{s,1}{s2,1}{ord,1};
                countavgL2 = countavgL2+1;
            elseif  strcmp(sublist_side{s},'R') && s==2
                ord = (no_cond+1)-i;
                ERPs{s,1}{s2,1}{ord,1} = double(squeeze(mean(EEG.data(:,:,find(events==cond(i))),3))); % list, subject, condition
                gavgR1 = gavgR1 + ERPs{s,1}{s2,1}{ord,1};
                countavgR1 = countavgR1+1;
            elseif strcmp(sublist_side{s},'R') && s==4
                ord = (no_cond+1)-i;
                ERPs{s,1}{s2,1}{ord,1} = double(squeeze(mean(EEG.data(:,:,find(events==cond(i))),3))); % list, subject, condition
                gavgR2 = gavgR2 + ERPs{s,1}{s2,1}{ord,1};
                countavgR2 = countavgR2+1;
            end
            %close all
            %plotdata = ERPs{s,1}{s2,1}{ord,1};
            %[maxval, maxidx] = max(abs(plotdata(:,:)),[],2);
            %[~, maxmaxidx] = max(maxval);
            %plottime = EEG.times(maxidx(maxmaxidx));
            %if plottime == EEG.times(end)
            %    plottime = EEG.times(end-1);
            %end
            %figure('Name',[subj '-' num2str(s)]);
            %timtopo(plotdata,chanlocs,...
            %    'limits',[EEG.times(1) EEG.times(end)],...
            %    'plottimes',plottime);
            %set(gcf,'Color','white');
            %pause
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

% plot grand average
for ii = 1:length(gavg);
    plotdata = gavg{ii,1};
    [maxval, maxidx] = max(abs(plotdata(:,:)),[],2);
    [~, maxmaxidx] = max(maxval);
    plottime = EEG.times(maxidx(maxmaxidx));
    if plottime == EEG.times(end)
        plottime = EEG.times(end-1);
    end
    figure('Name','grand average');
    timtopo(plotdata,chanlocs,...
        'limits',[EEG.times(1) EEG.times(end)],...
        'plottimes',plottime);
    set(gcf,'Color','white');
end

% identify limits for each peak
for ii = 1:length(gavg);
    gavgtmp = gavg{ii,1};
    %start_time = 0;
    %end_time = 200;

    for i = 1:dp
        mm_gavg(1,i) = max(gavgtmp(:,i)) - min(gavgtmp(:,i));
    end

    figure
    loc = zeros(1,no_peaks-1);
    peak = loc;
    for i=1:no_peaks-1
        peakdata = -mm_gavg(:,find(EEG.times==peaks_wind{i,1}):find(EEG.times==peaks_wind{i,2}));
        [peaks,locs] = max(peakdata);
        loc(1,i) = EEG.times(find(EEG.times==peaks_wind{i,1})-1+locs);
        peak(1,i) = peaks;
    end
    plot(EEG.times',mm_gavg);
    hold on;
    scatter(loc,-peak);

    if loc(1)~=0; loc = [0 loc]; end
    if loc(end)~=EEG.times(end);loc = [loc EEG.times(end)]; end
    for i = 1:length(loc)-1
        limit{i,1} = [loc(i) loc(i+1)];
    end
    limits{ii,1} = limit;
end

% identify top 5 electrodes for each peak from grand average, for each
% polarity


for ii = 1:length(gavg);
    limit = limits{ii,1};
    clear peaks locs
    peakdata = gavg{ii,1};
    
    for i = 1:no_peaks
        lim = limit{i,1};
        range = find(EEG.times==lim(1)):find(EEG.times==lim(2));
        e_peaks = zeros(el,1); 
        
        for e = 1:el
            [pospeak] = max(peakdata(e,range)); % greatest amplitude for each electrode, i.e. peak
            [negpeak] = min(peakdata(e,range)); % lowest amplitude for each electrode, i.e. trough
            if abs(pospeak) > abs(negpeak); e_peaks(e,1) = pospeak; else e_peaks(e,1) = negpeak;end
            e_peaks(e,2) = find(peakdata(e,:)==e_peaks(e,1));
        end

        [e_peaks_max ei] = sort(e_peaks(:,1),'descend'); % 
        e_locs_max = e_peaks(:,2);
        e_locs_max = e_locs_max(ei);
        e_max = e_peaks_max(1:no_ele);
        e_locs_max = e_locs_max(1:no_ele);
        [e_peaks_min ei] = sort(e_peaks(:,1),'ascend');
        e_locs_min = e_peaks(:,2);
        e_locs_min = e_locs_min(ei);
        e_min = e_peaks_min(1:no_ele);
        e_locs_min = e_locs_min(1:no_ele);

        e_sel=[];
        %if median(abs(e_max)) > median(abs(e_min))
            
            for j = 1:length(e_max)
                e_sel(1,j) = find(e_peaks(:,1)==e_max(j));
                ele{ii,1}{i,2}=median(e_locs_max);
                ele{ii,1}{i,4}=median(e_max);
            end
            ele{ii,1}{i,1}=e_sel;
            ele{ii,1}{i,3}=chanlocs(e_sel);
            
            for j = 1:length(e_min)
                e_sel(1,j) = find(e_peaks(:,1)==e_min(j));
                ele{ii,2}{i,2}=median(e_locs_min);
                ele{ii,2}{i,4}=median(e_min);
            end
            ele{ii,2}{i,1}=e_sel;
            ele{ii,2}{i,3}=chanlocs(e_sel);
            
        %else
            
        %    for j = 1:length(e_min)
        %        e_sel(1,j) = find(e_peaks(:,1)==e_min(j));
        %        ele{ii,1}{i,2}=median(e_locs_min);
        %        ele{ii,1}{i,4}=median(e_min);
        %    end
        %    ele{ii,1}{i,1}=e_sel;
        %    ele{ii,1}{i,3}=chanlocs(e_sel);
            
        %    for j = 1:length(e_max)
        %        e_sel(1,j) = find(e_peaks(:,1)==e_max(j));
        %        ele{ii,2}{i,2}=median(e_locs_max);
        %        ele{ii,2}{i,4}=median(e_max);
        %    end
        %    ele{ii,2}{i,1}=e_sel;
        %    ele{ii,2}{i,3}=chanlocs(e_sel);
        %end
        
        
    end
end


% plot electrodes for each peak
for ii = 1:length(gavg);
    peakdata = gavg{ii,1};
    figure
    subplot(4,7,1); topoplot([0 0 0 0 0], ele{ii,1}{1,3});
    subplot(4,7,2); topoplot([0 0 0 0 0], ele{ii,1}{2,3});
    subplot(4,7,3); topoplot([0 0 0 0 0], ele{ii,1}{3,3});
    subplot(4,7,4); topoplot([0 0 0 0 0], ele{ii,1}{4,3});
    subplot(4,7,5); topoplot([0 0 0 0 0], ele{ii,1}{5,3});
    subplot(4,7,6); topoplot([0 0 0 0 0], ele{ii,1}{4,3});
    subplot(4,7,7); topoplot([0 0 0 0 0], ele{ii,1}{5,3});
    subplot(4,7,8); topoplot(peakdata(:,ele{ii,1}{1,2}), chanlocs,'maplimits',[-1 1]);
    subplot(4,7,9); topoplot(peakdata(:,ele{ii,1}{2,2}), chanlocs,'maplimits',[-5 5]);
    subplot(4,7,10); topoplot(peakdata(:,ele{ii,1}{3,2}), chanlocs,'maplimits',[-5 5]);
    subplot(4,7,11); topoplot(peakdata(:,ele{ii,1}{4,2}), chanlocs,'maplimits',[-3 3]);
    subplot(4,7,12); topoplot(peakdata(:,ele{ii,1}{5,2}), chanlocs,'maplimits',[-7 7]);
    subplot(4,7,13); topoplot(peakdata(:,ele{ii,1}{4,2}), chanlocs,'maplimits',[-4 4]);
    subplot(4,7,14); topoplot(peakdata(:,ele{ii,1}{5,2}), chanlocs,'maplimits',[-7 7]);
    subplot(4,7,15); topoplot([0 0 0 0 0], ele{ii,2}{1,3});
    subplot(4,7,16); topoplot([0 0 0 0 0], ele{ii,2}{2,3});
    subplot(4,7,17); topoplot([0 0 0 0 0], ele{ii,2}{3,3});
    subplot(4,7,18); topoplot([0 0 0 0 0], ele{ii,2}{4,3});
    subplot(4,7,19); topoplot([0 0 0 0 0], ele{ii,2}{5,3});
    subplot(4,7,20); topoplot([0 0 0 0 0], ele{ii,2}{4,3});
    subplot(4,7,21); topoplot([0 0 0 0 0], ele{ii,2}{5,3});
    subplot(4,7,22); topoplot(peakdata(:,ele{ii,2}{1,2}), chanlocs,'maplimits',[-1 1]);
    subplot(4,7,23); topoplot(peakdata(:,ele{ii,2}{2,2}), chanlocs,'maplimits',[-5 5]);
    subplot(4,7,24); topoplot(peakdata(:,ele{ii,2}{3,2}), chanlocs,'maplimits',[-5 5]);
    subplot(4,7,25); topoplot(peakdata(:,ele{ii,2}{4,2}), chanlocs,'maplimits',[-3 3]);
    subplot(4,7,26); topoplot(peakdata(:,ele{ii,2}{5,2}), chanlocs,'maplimits',[-7 7]);
    subplot(4,7,27); topoplot(peakdata(:,ele{ii,2}{4,2}), chanlocs,'maplimits',[-4 4]);
    subplot(4,7,28); topoplot(peakdata(:,ele{ii,2}{5,2}), chanlocs,'maplimits',[-7 7]);
end

% topoplots for specific latency window
%lat = [144]; %ms
%lat = [272:4:324]; %ms
%lat = [172:4:292]; %ms
%lat = [56]; %ms
lat = [128:4:136]; %ms
lat = (lat+200)/4;
    figure
for ii = 1:length(gavg);
    peakdata = gavg{ii,1};
    subplot(4,1,ii); topoplot(mean(peakdata(:,lat),2), chanlocs,'maplimits',[-5 5]);
end

peakselect = 4;
elegrpselect = 1;
figure
col = {'b','b--','r','r--'};
same_elegrp = 2;
times = find(((EEG.times >= -100).*(EEG.times <= 500))==1);
for ii = [1 3]; 
    peakdata = gavg{ii,1};
    if same_elegrp>0
        ERP = mean(peakdata(ele{same_elegrp,elegrpselect}{peakselect,1},:),1);
    else
        ERP = mean(peakdata(ele{ii,elegrpselect}{peakselect,1},:),1);
    end
    plot(EEG.times(times),ERP(times),col{ii}); hold on;
end
%legend('Healthy controls', 'healthy unaff', 'CRPS patients','patient unaff')
legend('Healthy controls', 'CRPS patients','patient unaff');
ylabel('Amplitude, uV') % label for y axis
xlabel('Time (ms)') % label for x axis
set(gca,'FontSize',15);
set(findall(gcf,'type','text'),'FontSize',15);

% identify peaks for each subject and condition
LAT = NaN('double');
AMP = NaN('double');

for s = 1:length(subjects)
    for s2 = 1:length(subjects{s,1}) 
        affs = aff_side(s2);
        for i = 1:length(ERPs{s,1}{s2,1}) % for 5 conditions
            ERP = ERPs{s,1}{s2,1}{i,1};
            %figure
            %peakdata = ERP;
            %plot(EEG.times',peakdata);hold on;
            for p = 1:length(limits{1,1}) % peak
                for iv=1:2 % elegrp hand represented by ii at 1 (left) or 2 (right)
                    if s==1 || s==2
                        if iv==1; ii=1;
                        else ii = 2;
                        end
                    else
                        if iv==1; ii=3;
                        else ii = 4;
                        end
                    end
                    %if strcmp(sublist_side{s},'L'); ii=1; elseif strcmp(sublist_side{s},'R') ii=2; end % which side is it?
                    lim = limits{ii,1}{p,1};
                    range = find(EEG.times==lim(1)):find(EEG.times==lim(2)); 
                    for iii = 1:2 % elegrp
                        ERPer = mean(ERP(ele{ii,iii}{p,1},range),1);
                        GFPer = std(ERP(:,range),1);
                        if ~isempty(ele_left) && ~isempty(ele_right)
                            if iv==1 && affs==1; ele_sel = ele_left;
                            elseif iv==1 && affs==2; ele_sel = ele_right;
                            elseif iv==2 && affs==1; ele_sel = ele_right;
                            elseif iv==2 && affs==2; ele_sel = ele_left;
                            end
                            ERPreg = mean(ERP(ele_sel,range),1);
                        else ERPreg = ERPer;
                        end
                        if ele{ii,iii}{p,4}>0
                            [pospeak locs] = max(ERPer);
                            if pospeak<0; [pospeak locs] = min(ERPer);end;
                            actpeak = pospeak; %list, subject, condition, peak, EleGroup
                        else
                            [negpeak locs] = min(ERPer);
                            if negpeak>0; [negpeak locs] = max(ERPer);end;
                            actpeak = negpeak;
                        end
                        
                        [gfppeak locs] = max(GFPer);
                        LAT(s,s2,i,p,iv,iii) = EEG.times(range(1)-1+find(ERPer==actpeak));  
                        %AMP(s,s2,i,p,iv,iii) = squeeze(mean(ERPer(1,max(1,(find(ERPer==actpeak)-floor(((lim(2)-lim(1))/4)/8))):min(length(ERPer),(find(ERPer==actpeak)+floor(((lim(2)-lim(1))/4)/8)))),2)); % range divided by 8
                        %AMP(s,s2,i,p,iv,iii) = squeeze(mean(ERPer(1,max(1,(find(GFPer==gfppeak)-floor(12/4))):min(length(ERPer),(find(GFPer==gfppeak)+floor(12/4)))),2)); % range divided by 8
                        AMP(s,s2,i,p,iv,iii) = squeeze(mean(ERPer(1,:),2)); % range divided by 8
                        GFP(s,s2,i,p) = squeeze(mean(GFPer(1,max(1,(find(GFPer==gfppeak)-floor(12/4))):min(length(GFPer),(find(GFPer==gfppeak)+floor(12/4)))),2)); % range divided by 8
                        
                        %plot(EEG.times, mean(ERP(ele{ii,iii}{p,1},:),1)); 
                        %hold on;
                        %scatter(LAT(s,s2,i,p,iv,iii),AMP(s,s2,i,p,iv,iii));
                        %pause;
                        %close all;
                        
                        %topoplot(ERP(:,find(EEG.times==LAT(s,s2,i,p,iv,iii))), chanlocs);
                        %title(['hand ' num2str(iv) '; peak ' num2str(p)]);
                        %figure
                        %topoplot([0 0 0 0 0], ele{iv,iii}{p,3});
                        %datap = find(EEG.times==LAT(s,s2,i,p,iv,iii));
                        %scatter(LAT(s,s2,i,p,iv,iii)*ones(size(ERP,1),1),ERP(:,datap));
                    end
                    %hold off;             
            
                end
            end
            
        end
        %topoplot(ERP(:,find(EEG.times==LAT(s,s2,i,p,iv,iii))), chanlocs);
        %title(['hand ' num2str(iv) '; peak ' num2str(p)]);
        %pause
        %close all;
    end
end

figure; 
ss= [2 4];
for s2 = 1:length(subjects{s,1}) 
    ERP_right = zeros(92,250);
    ERP_left = zeros(92,250);
    countright = 0;
    countleft = 0;
    for s =ss
        
        
        for i = 1:length(ERPs{s,1}{s2,1})
            ERP = ERPs{s,1}{s2,1}{i,1};
            for p = 4 % peak
                
                for iv=1:2 % hand
                    if strcmp(sublist_side{s},'L'); 
                        ii=1; 
                        lim = limits{ii,1}{p,1};
                        range_left = find(EEG.times==lim(1)):find(EEG.times==lim(2)); 
                        ERP_left = ERP_left+ERP;
                        countleft=countleft+1;
                    elseif strcmp(sublist_side{s},'R') 
                        ii=2; 
                        lim = limits{ii,1}{p,1};
                        range_right = find(EEG.times==lim(1)):find(EEG.times==lim(2)); 
                        ERP_right = ERP_right+ERP;
                        countright=countright+1;
                    end % which side is it?
                    
                end
            end
            
        end
       
    end
    ERP_left = ERP_left/countleft;
    ERP_right = ERP_right/countright;
    ind = s2*2-1;
    subplot(13,2,ind); topoplot(mean(ERP_left(:,range_left),2), chanlocs);
    subplot(13,2,ind+1); topoplot(mean(ERP_right(:,range_right),2), chanlocs);
end

% create results tables for export to excel / SPSS
usublist = unique(sublist','stable');
nrow = length(usublist);
AMP2 = reshape(AMP,size(AMP,1),size(AMP,2),no_cond,no_peaks,no_elegrp);
LAT2 = reshape(LAT,size(LAT,1),size(LAT,2),no_cond,no_peaks,no_elegrp);
ncol = 2*size(AMP2,3)*size(AMP2,4)*size(AMP2,5); % hand, condition, peak, EleGroup
nrowhead = 1;
ncolhead = 1;
amp_results = cell(nrow+ncolhead,ncol+nrowhead);
amp_results(ncolhead+1:end,1) = usublist;
for p = 1:no_peaks
    for e = 1:no_elegrp
        for c = 1:no_cond
            for h = 1:2
                ind = ((p-1)*no_cond*2*no_elegrp + (e-1)*2*no_cond + (h-1)*no_cond + c);
                amp_results{1,nrowhead+ind} = [peaks_nme{p} '_' elegrp_nme{e} '_' hand_nme{h} '_' cond_nme{c}];
            end
        end
     end
end

lat_results = amp_results;

AMP3 = reshape(AMP2,size(AMP2,1)/2,size(AMP2,1)/2,size(AMP2,2),no_cond,no_peaks,no_elegrp); % hand, group, subject, cond, peaks, elegrp
AMP4 = permute(AMP3,[3 2 4 1 6 5]); % subject, group, cond, hand, elegrp, peaks 
AMP5 = reshape(AMP4,size(AMP4,1)*size(AMP4,2),2*no_cond*no_peaks*no_elegrp);
amp_results(ncolhead+1:end,nrowhead+1:end) = num2cell(AMP5);

LAT3 = reshape(LAT2,size(LAT2,1)/2,size(LAT2,1)/2,size(LAT2,2),no_cond,no_peaks,no_elegrp);
LAT4 = permute(LAT3,[3 2 6 5 4 1]);
LAT5 = reshape(LAT4,size(LAT4,1)*size(LAT4,2),2*no_cond*no_peaks*no_elegrp);
lat_results(ncolhead+1:end,nrowhead+1:end) = num2cell(LAT5);

% create latency table for SPM source analysis

latencies = squeeze(mean(LAT4,3)); %mean elegrp
latencies = squeeze(mean(latencies,4)); %mean cond
latencies = permute(latencies,[1 4 2 3]); % subject, group, peaks, hand
latencies = reshape(latencies,size(latencies,1)*size(latencies,2)*2,no_peaks); % group * hand * subject, condition * peak
save latencies_unflipped.dat latencies -ascii

% create matrices with meaned values

%AMP_allfingers = squeeze(mean(AMP4,5));
%AMP_allfingers_peak1 = squeeze(AMP_allfingers(:,:,:,1,:));
%AMP_allfingers_peak2 = squeeze(AMP_allfingers(:,:,:,2,:));
%AMP_allfingers_peak3 = squeeze(AMP_allfingers(:,:,:,3,:));
%AMP_allfingers_peak4 = squeeze(AMP_allfingers(:,:,:,4,:));%
%AMP_allfingers_peak5 = squeeze(AMP_allfingers(:,:,:,5,:));
%AMP_allfingers_peak1 = reshape(AMP_allfingers_peak1,size(AMP4,1)*size(AMP4,2),2*no_elegrp);
%AMP_allfingers_peak2 = reshape(AMP_allfingers_peak2,size(AMP4,1)*size(AMP4,2),2*no_elegrp);
%AMP_allfingers_peak3 = reshape(AMP_allfingers_peak3,size(AMP4,1)*size(AMP4,2),2*no_elegrp);
%AMP_allfingers_peak4 = reshape(AMP_allfingers_peak4,size(AMP4,1)*size(AMP4,2),2*no_elegrp);
%AMP_allfingers_peak5 = reshape(AMP_allfingers_peak5,size(AMP4,1)*size(AMP4,2),2*no_elegrp);
%AMP_allfingers_peaks = cat(2,AMP_allfingers_peak1,AMP_allfingers_peak2,AMP_allfingers_peak3,AMP_allfingers_peak4,AMP_allfingers_peak5);
%AMP_peak_group = squeeze(mean(AMP4,3));
%AMP_peak_group = squeeze(mean(AMP_peak_group,4));
%AMP_peak_group = squeeze(mean(AMP_peak_group,4));

AMP_peaks_elegrp = squeeze(mean(AMP4,3));
AMP_peaks_elegrp = squeeze(mean(AMP_peaks_elegrp,3));
AMP_peaks_elegrp = reshape(AMP_peaks_elegrp,size(AMP4,1)*size(AMP4,2),no_peaks*no_elegrp);

AMP_peaks_elegrp_hand = squeeze(mean(AMP4,3));
AMP_peaks_elegrp_hand = reshape(AMP_peaks_elegrp_hand,size(AMP4,1)*size(AMP4,2),size(AMP4,4)*no_peaks*no_elegrp);
ncol = size(AMP_peaks_elegrp_hand,2); % hand, condition, peak, EleGroup
nrowhead = 1;
ncolhead = 1;
amp_results_nocond = cell(nrow+ncolhead,ncol+nrowhead);
amp_results_nocond(ncolhead+1:end,1) = usublist;
for p = 1:no_peaks
    for e = 1:no_elegrp
        %for c = 1:no_cond
            for h = 1:2
                ind = ((p-1)*2*no_elegrp + (e-1)*2 + h);
                amp_results_nocond{1,nrowhead+ind} = [peaks_nme{p} '_' elegrp_nme{e} '_' hand_nme{h}];
            end
        %end
     end
end
amp_results_nocond(ncolhead+1:end,nrowhead+1:end) = num2cell(AMP_peaks_elegrp_hand);
