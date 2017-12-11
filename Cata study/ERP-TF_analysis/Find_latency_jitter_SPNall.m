%---select data to load---%
%select = 'TF'; TFmethod = '-FT'; % or, '-EL'
close all
clear all
Dpath = 'C:\Data\Catastrophising study\Preprocessed';
cd(Dpath)
select = 'ERP'; TFmethod = '';
%savenametype='';
%for i = use_etype
    savenametype = '_1_3_5_7_a_a_a_a_b_b_b_b';
%end
handtype = {'R'}; hand_ana=1; % 1:LR, 2:AU
hand_nme={'right'};
load([select TFmethod savenametype '_data_SPN.mat'])
eventtypes = {'c1','c3','c5','c7','c2a','c4a','c6a','c8a','c2b','c4b','c6b','c8b'}; use_etype = [1 2 3 4; 5 6 7 8; 9 10 11 12]; no_cond = size(use_etype);
benchmark = 3;
modify=[1 2];

t = [-2500,-2000]; %ms
dp = (t+3000)/2;
baseline = [1 250];

anatype = {'evoked', 'induced','itc'}; anaDAT=1; % 1:evoked, 2:induced, 3: ITC
min_trials_gavg =10; % minimum no. of trials per condition
output_ntrials=0;
data_max=0;
%baseline = [1400 1500];
%baseline = [2500 2750];
%--------------%

gavg = cell(size(use_etype,1),2);
for c = 1:size(use_etype,1)
    gavg{c,1} = zeros(el,length(times),nf);
    gavg{c,2} = 0;
    for g = 1:length(subjects)
        for s = 1:length(subjects{g,1}) 
            for i = use_etype(c,:)
                if DATs{g,1}{s,1}{i,6}>=min_trials_gavg
                    gavg{c,1} = gavg{c,1} + DATs{g,1}{s,1}{i,anaDAT};
                    gavg{c,2} = gavg{c,2}+1;
                end
            end
        end
    end
    gavg{c,1} = gavg{c,1}/gavg{c,2};
end

r=[];
rm=[];
lm = length(modify);
for ng = 1:lm
    for e = 1:size(gavg{modify(ng),1},1)
        dat1 = gavg{modify(ng),1}(e,dp(1):dp(2))';
        for d = 1:(dp(2)-dp(1))
            datSUB = gavg{benchmark,1}(e,(dp(1):dp(2))-d+1)';
            datADD = gavg{benchmark,1}(e,(dp(1):dp(2))+d)';
            r(ng,e,d,1) = corr(dat1,datSUB,'type','Spearman');
            r(ng,e,d,2) = corr(dat1,datADD,'type','Spearman');
        end
    end
    rm(:,ng,:) = squeeze(mean(r(ng,:,:,:),2));
    [rmaxSUB,imaxSUB] = max(rm(:,ng,1));
    [rmaxADD,imaxADD] = max(rm(:,ng,2));
    if rmaxSUB>rmaxADD
        ms_jitter(ng) = imaxSUB*2
    else
        ms_jitter(ng) = -imaxADD*2
    end
end

figure;plot(gavg{benchmark,1}(15,:),'k');hold on;plot(gavg{modify(1),1}(15,:),'r');hold on;plot(gavg{modify(2),1}(15,:),'r--')