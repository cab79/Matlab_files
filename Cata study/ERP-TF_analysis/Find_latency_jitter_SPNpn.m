%---select data to load---%
%select = 'TF'; TFmethod = '-FT'; % or, '-EL'
close all
clear all
Dpath = 'C:\Data\Catastrophising study\Preprocessed';
cd(Dpath)
select = 'ERP'; TFmethod = '';
%savenametype='';
%for i = use_etype
    savenametype = '_a_b_a_b_a_b_a_b';
%end
handtype = {'R'}; hand_ana=1; % 1:LR, 2:AU
hand_nme={'right'};
load([select TFmethod savenametype '_data_SPNpn.mat'])
eventtypes = {'c2a','c2b','c4a','c4b','c6a','c6b','c8a','c8b'}; use_etype = [1 3 5 7; 2 4 6 8]; no_cond = size(use_etype);

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
for e = 1:size(gavg{2,1},1)
    dat1 = gavg{2,1}(e,dp(1):dp(2))';
    for d = 1:(dp(2)-dp(1))
        dat2 = gavg{1,1}(e,(dp(1):dp(2))-d+1)';
        r(e,d) = corr(dat1,dat2,'type','Spearman');
    end
end
r = mean(r,1);
[rmax,imax] = max(r);
ms_jitter = imax*2