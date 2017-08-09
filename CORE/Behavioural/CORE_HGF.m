clear all
dbstop if error
addpath('C:\Data\Matlab\Matlab_files\CORE\Experiment');
addpath('C:\Data\Matlab\HGF\HGFv5.0');
addpath('C:\Data\Matlab\Matlab_files\CORE\Behavioural\HGF_functions')

dname='C:\Data\CORE\Design';
rname='C:\Data\CORE\Behaviour';
aname='C:\Data\CORE\Behaviour\July2017\HGF_Results\SDTnp_mm_thresh_b0_sd05_0-150'; bayes_opt=0; part='part4';
%aname='C:\Data\CORE\Behaviour\July2017\HGF_Results\KF_mmpredict_soft'; bayes_opt=0; part='part4';
%aname='C:\Data\CORE\Behaviour\July2017\HGF_Results\3lev_soft_a2'; bayes_opt=0; part='part2';
%aname='C:\Data\CORE\Behaviour\July2017\HGF_Results\SDT_Bayes_part4'; bayes_opt=0; part='part2';
%cd(rname);

% directory containing cosmo projections - part4 analysis only
cosdir = 'C:\Data\CORE\cosmo';
% range of cosdir extensions, each with different MM projections
cosdirext = {'LDA_part4_timechan_0_150'};
% generic filename
cosname = '_4_mm_projection.mat';
% options
mm_trials = 0; % 1 = use only mismatch trials
mm_positive = 0; % 1 = use only positive values of mm projection

overwrite=0;
hgf=1;
sdt=0;
acc = 0;

% input codes (for stimulus types 1 to 7)
% code 'NaN' to ignore, otherwise give it a unique value
%ic = [1 2 3 4 5 6 7]; prc_config = 'tapas_hgf_binary_pu_CORE_config'; % "alpha 7"
%ic = [1 2 3 4 NaN NaN NaN]; prc_config = 'tapas_hgf_binary_pu_CORE_config_alpha4'; % "alpha 4"
%ic = [1 1 2 2 NaN NaN NaN]; prc_config = 'tapas_hgf_binary_pu_CORE_config_alpha2';% "alpha 2a"
%ic = [1 2 1 2 NaN NaN NaN]; prc_config = 'GBM_config_3lev'; obs_config = 'logrt_softmax_binary_soft_config'; nstim=[];% either muin or alphain set to 4
%ic = [1 1 1 1 NaN NaN NaN]; prc_config = 'tapas_hgf_binary_pu_CORE_config_alpha1';% 
%ic = [1 1 1 1 NaN NaN NaN]; prc_config = 'tapas_sdt_config';% 
%ic = [1 2 3 4 NaN NaN NaN]; prc_config = 'GBM_config_SDT_noprior'; obs_config = 'logrt_softmax_binary_RTsoft_config'; nstim=[];% either muin or alphain set to 4
%ic = [1 1 1 1 NaN NaN NaN]; prc_config = 'GBM_config_3lev'; obs_config = 'tapas_bayes_optimal_binary_config'; nstim=[];
%ic = [1 1 1 1 NaN NaN NaN]; prc_config = 'GBM_config_SDT_noprior'; obs_config = 'logrt_softmax_binary_soft_config'; nstim=[];
ic = [1 2 1 2 NaN NaN NaN]; prc_config = 'GBM_config_SDT_noprior'; obs_config = 'mismatch_response_binary_sd05_config'; nstim=[];
%ic = [1 1 1 1 NaN NaN NaN]; prc_config = 'GBM_config_KF'; obs_config = 'mismatch_softmax_binary_config'; nstim=[];

% load .xlsx file containing 'Participant_ID', 'Group', and covariates
pdatfile = 'C:\Data\CORE\Participant_data.xlsx';
% names of headers in the above xls file:
    subhead = 'Subject';
    grphead = 'Group';
    inchead = 'Include_behav';
    CRPSsidehead = 'CRPSlr';
% which codes to analyse in 'Include' columns in participant data file?
include_codes = [1];

% event numbers for each condition (1st = change, 2nd = no change)
cond1 = [1 3]; %0.1 prob, 1 digit change, left hand
cond2 = [2 4]; %0.1 prob, 3 digit change, left hand
cond3 = [5 7]; %0.1 prob, 1 digit change, right hand
cond4 = [6 8]; %0.1 prob, 3 digit change, right hand
cond5 = [9 11]; %0.3 prob, 1 digit change, left hand
cond6 = [10 12]; %0.3 prob, 3 digit change, left hand
cond7 = [13 15]; %0.3 prob, 1 digit change, right hand
cond8 = [14 16]; %0.3 prob, 3 digit change, right hand
cond9 = [17 19]; %0.5 prob, 1 digit change, left hand
cond10 = [18 20]; %0.5 prob, 3 digit change, left hand
cond11 = [21 23]; %0.5 prob, 1 digit change, right hand
cond12 = [22 24]; %0.5 prob, 3 digit change, right hand

factors = {'CP', 'Side', 'DC'}; % oly need to label those that 
levels = {{'10' '30' '50'}, {'L' 'R'}, {'1' '3'}}; % oly need to label those that 
levels_aff = {'aff' 'unaff'}; analyse_aff=1;
factor_matrix = [
          1 1 1
          1 1 2
          1 2 1
          1 2 2
          2 1 1
          2 1 2
          2 2 1
          2 2 2
          3 1 1
          3 1 2
          3 2 1
          3 2 2
          ];

%% RUN

if ~exist(aname,'dir') && hgf
    mkdir(aname)
end

if analyse_aff
    levels{2} = levels_aff;
    opt_name = '_aff';
else
    opt_name = '';
end

header={'Subject'};
for fa = 1:length(factor_matrix)
    header = horzcat(header,[factors{1} levels{1}{factor_matrix(fa,1)} '_' factors{2} levels{2}{factor_matrix(fa,2)} '_' factors{3} levels{3}{factor_matrix(fa,3)}]);
end
acc_all=header;
rt_all=header;
dp_all=header;
cri_all=header;
hit_all=header;
fa_all=header;
dp_all(1,end+1)= {'total'};
cri_all(1,end+1)= {'total'};
hit_all(1,end+1)= {'total'};
fa_all(1,end+1)= {'total'};

[~,~,pdata] = xlsread(pdatfile);
grp_col = find(strcmp(pdata(1,:),grphead));
sub_col = find(strcmp(pdata(1,:),subhead));
inc_col = find(strcmp(pdata(1,:),inchead));
side_col = find(strcmp(pdata(1,:),CRPSsidehead));
inc_idx = cellfun(@(x) ismember(x,include_codes), pdata(2:end,inc_col), 'UniformOutput', 0);
inc_idx = find(cell2mat(inc_idx));
CRPSsides = pdata(2:end,side_col);

% add random sides to the healthy subjects (who have no CRPS side)
sidenan = cell2mat(cellfun(@(x) isnan(x), CRPSsides, 'UniformOutput', 0));
sidepool = CRPSsides(~sidenan);
sn_idx=find(sidenan);
for s=1:length(sn_idx)
    CRPSsides(sn_idx(s)) = sidepool(randi([1 numel(sidepool)]));
end

files_ana = inc_idx;
for f = files_ana'
    
    C = strsplit(pdata{f+1,sub_col},'CORE');
    
    file = dir(fullfile(dname,['dt_' C{2} '_' part '*']));
    if isempty(file);continue;end;
    dt_name = file.name
    load(fullfile(dname,dt_name)); 
    design=dt.design;
    u = design(3,:)'; % change = 1, no change = 0;
    comb = [u'];
    
    if strcmp(part,'part2') && (acc || sdt)
        acc_all(1+f,1)= pdata(f+1,sub_col);
        rt_all(1+f,1)= pdata(f+1,sub_col);
        dp_all(1+f,1)= pdata(f+1,sub_col);
        cri_all(1+f,1)= pdata(f+1,sub_col);
        hit_all(1+f,1)= pdata(f+1,sub_col);
        fa_all(1+f,1)= pdata(f+1,sub_col);
    end
        
    if strcmp(part,'part2')
        RT_name = ['RT_' dt_name(4:end)];
        load(fullfile(rname,RT_name));
        comb = [comb;RT(1,:)];
    end
    
    [hand,dc,cp,bi,blockii,btypes] = blocktype(dname,dt_name);
    
    % swap sides for people who are right affected
    if analyse_aff && strcmp(CRPSsides{f},'R')
        i1 = hand==1;
        i2 = hand==2;
        hand(i1)=2;
        hand(i2)=1;
        temp = cell2mat(btypes(:,2));
        i1 = temp==1;
        i2 = temp==2;
        temp(i1)=2;
        temp(i2)=1;
        btypes(:,2) = num2cell(temp)';
        i1 = bi==1;
        i2 = bi==2;
        i3 = bi==3;
        i4 = bi==4;
        i5 = bi==5;
        i6 = bi==6;
        i7 = bi==7;
        i8 = bi==8;
        i9 = bi==9;
        i10 = bi==10;
        i11 = bi==11;
        i12 = bi==12;
        bi(i1)=3;
        bi(i2)=4;
        bi(i3)=1;
        bi(i4)=2;
        bi(i5)=7;
        bi(i6)=8;
        bi(i7)=5;
        bi(i8)=6;
        bi(i9)=11;
        bi(i10)=12;
        bi(i11)=9;
        bi(i12)=10;
    end

    
    for i = 2:size(design,2) % ignore first trial as unlikely to incur a response
        if design(2,i)==0 % mark changes of block as being a stimulus change
            comb(1,i)=1;
        end
    end
    
    % create labels for block transitions between hands and within hands
    u2=zeros(1,size(comb,2));
    for i = 2:length(blockii) 
        if hand(i) ~= hand(i-1) % label transitions between hands as 7
            u2(blockii(i))=ic(7); 
        elseif dc(i) ~= dc(i-1) % label transitions within hands (from DC1 to DC3) as 5 for left, 6 for right
            if hand(i)==1
                u2(blockii(i))=ic(5); 
            elseif hand(i)==2
                u2(blockii(i))=ic(6);
            end
        end
    end
    comb(1,blockii)=u2(1,blockii);
    for i = 1:length(u2)
        if design(2,i)~=0
            ii = find(any(cell2mat(btypes(:,1))==design(2,i),2));
            handval = btypes{ii,2};
            dcval = btypes{ii,3};
            if handval==1 && dcval==1
                u2(i)=ic(1);
            elseif handval==1 && dcval==2
                u2(i)=ic(2);
            elseif handval==2 && dcval==1
                u2(i)=ic(3);
            elseif handval==2 && dcval==2
                u2(i)=ic(4);
            end
        end
    end
    
    % CHANGE REMAINING ZEROS IN U2 to be same as previous trial
    zi = find(u2==0);
    zi(zi==1)=[]; % ignore first trial
    u2(zi) = u2(zi-1);
    u2(1) = u2(2);
    
    u=[];
    u(1,:) = comb(1,:)>0;
    %u2(u2>2) =1;
    u(2,:) = u2;
    u=u';
    u(isnan(u(:,2)),1)=NaN;
    
    transi = ismember(u2,[1:4]); %use only transitions of interest (1-4 from u2)
    % create index of cond numbers corresponding to each trial
    condi =nan(1,length(u2));
    blockii_end = [blockii length(u2)+1];
    for b = 1:length(bi)
        condi(blockii_end(b):blockii_end(b+1)-1) = bi(b);
    end
    conds = sort(unique(bi));
    
    % Response data
    if strcmp(part,'part2')
        ISI = 1.0;
        thresh = 0.2; % min RT that is realistic
        num_trial_lag = [1 1 1]; % number of trials over which response is allowed to lag for each of 0.1 prob, 0.3 prob, 0.5 prob.

        for i = 1:size(design,2)-1
            if any(cond1 == design(2,i)) || any(cond2 == design(2,i)) || any(cond3 == design(2,i)) || any(cond4 == design(2,i)) 
                lag = num_trial_lag(1);
            elseif any(cond5 == design(2,i)) || any(cond6 == design(2,i)) || any(cond7 == design(2,i)) || any(cond8 == design(2,i)) 
                lag = num_trial_lag(2);
            elseif any(cond9 == design(2,i)) || any(cond10 == design(2,i)) || any(cond11 == design(2,i)) || any(cond12 == design(2,i)) 
                lag = num_trial_lag(3);
            end

            comb(3,i)=0; % row 3 is the corrected response time on the corrected trial
            if comb(1,i)>0 % if the stimulus actually changed, but the response occured on the next trial, count it as occuring on the current trial but with a longer RT
                for j = 1:lag+1
                    move=0;
                    if comb(3,i)==0 && comb(2,i+(j-1))>0 && comb(3,i-1)~=ISI*1+comb(2,i+(j-1)) && comb(2,i+(j-1))+(j-1)*ISI>thresh
                        if comb(1,i+(j-1))==0 || j==1 % if we are considering the current trial only, or if there is no stimulus on the next trial, then it's always ok to register the response on te current trial.
                            move=1;
                        elseif comb(1,i+(j-1))>0 && j>1 %if there is a stimulus on the next trial AND
                            if comb(2,i+(j-1))>0 && comb(2,i+(j-1))<thresh % the response is less than threshold OR
                                move=1;
                            elseif (comb(1,i+j)==0 && comb(2,i+j)>0) % there is no stim on the subsequent trial but a response
                                move=1;
                            end
                        end
                    end
                    if move==1
                        comb(3,i)=comb(2,i+(j-1)) + (j-1)*ISI; % adds on the ISI from RT on next trial and places it on the current trial
                    end
                end
            end
            if comb(1,i)==0 && comb(2,i)>0 && i>1 % if no stimulus change, but a response, 
                if ~any(comb(1,i-(j-1):i)>0) % and there was no stim change on previous trial,, so it can't be a delayed response to a stim change
                    if comb(2,i)>thresh;
                        comb(3,i)=comb(2,i); % count it as a false positive response
                    elseif comb(3,i-1)==0
                        comb(3,i-1)=comb(2,i)+1*ISI; % or a response on the previous trial if too fast for current trial
                    end
                end
            end
        end 
        y = double(comb(3,:)>0); % BINARY response
        y(2,:) = comb(3,:); % RTs
        rtnan = y(2,:);
        rtnan(rtnan==0)=nan;
        y(2,:) = rtnan;
        y=y';
        
        %% analysis % correct and RTs for each condition 1-12
        pc_corr = nan(length(conds),1);
        rt_corr = nan(length(conds),1);
        for c = conds
            uyr=[];
            uyr(1,:) = u(condi==c,1);
            uyr(2,:) = y(condi==c);
            uyr(3,:) = comb(3,condi==c);
            uyr(4,:) = transi(condi==c);
            uyr = uyr(:,uyr(4,:)==1);
            pc_corr(c) = 100*sum(double(uyr(1,:)==uyr(2,:)))/size(uyr,2);
            hit(c) = length(intersect(find(uyr(1,:)==1),find(uyr(2,:)==1)))/length(find(uyr(1,:)==1)); % hit rate for SDT
            fa(c) = length(intersect(find(uyr(1,:)==0),find(uyr(2,:)==1)))/length(find(uyr(1,:)==0)); % false alarm rate for SDT
            rt_corr(c) = mean(uyr(3,intersect(find(uyr(1,:)==1), find(uyr(2,:)==1))));
        end

        % create data spreadsheets
        acc_all(f+1,2:size(pc_corr,1)+1)=num2cell(pc_corr');
        rt_all(f+1,2:size(rt_corr,1)+1)=num2cell(rt_corr');

        if sdt
            hit(c+1) = length(intersect(find(u(:,1)==1),find(y==1)))/length(find(u(:,1)==1)); % hit rate for SDT
            fa(c+1) = length(intersect(find(u(:,1)==0),find(y==1)))/length(find(u(:,1)==0)); % false alarm rate for SDT
            [dp,cri] = signal_detection(hit,fa);
            dp_all(f+1,2:size(dp,2)+1)=num2cell(dp');
            cri_all(f+1,2:size(cri,2)+1)=num2cell(cri');
            hit_all(f+1,2:size(dp,2)+1)=num2cell(hit');
            fa_all(f+1,2:size(cri,2)+1)=num2cell(fa');
        end

    elseif strcmp(part,'part4')
        load(fullfile(cosdir,cosdirext{1},['CORE' C{2} cosname]));
        y=nan(length(u),3);
        
        %mm=class_proj(1,:)';
        
        % standardise mm
        %mm=zscore(mm);
        
        % add raw mm to column 3 of y
        y(tnums,3)=mm;
        if mm_trials
            y(:,3) = y(:,3).*u(:,1);
        end
        if mm_positive
            y(:,3) = y(:,3).*(y(:,3)>0);
        end
        y(y(:,3)==0,3) = nan;
        
        % class predictions for softmax function
        %mm_thresh = 0; % number of SDs
        %class1=1*single(mm>mm_thresh);%CAB 
        %class2=0*single(mm<=mm_thresh);%CAB 
        %predicted=class1+class2;%CAB 
        %y(tnums,4)=predicted;
    end
    

    %conds_ana = [1 3 5];
    %conds_ana = [2 4 6];
    %idx=[];
    %for i = 1:length(conds_ana)
    %    eval(['conds = cond' num2str(conds_ana(i)) ';']); 
    %    for j = 1:length(conds)
    %        idx = [idx find(design(2,:)==conds(j))];
    %    end
    %end
    %idx = sort(idx);
    %design = design(:,idx);
    %comb = comb(:,idx); 
    

    %for i = 1:length(u)
    %    if u(i)==0.1
    %        u(i)=1.75;
    %    end
    %end
    
    

    %% HGF
    % prc: perceptual; obs:observation; opt:optimisation
    prc_model = prc_config;
    %prc_model = 'tapas_hgf_binary_pu_config';
    obs_model = obs_config;
    %obs_model = 'tapas_bayes_optimal_binary_config'; %BAYES OPTIMAL
    opt_algo = 'tapas_quasinewton_optim_config';
    
    sname = [C{2} '_bopars' opt_name '.mat'];
    if hgf && exist(fullfile(aname,sname),'file') && overwrite==0
        load(fullfile(aname,sname));
    elseif hgf
        if isempty(nstim)
            nst=length(u);
        else
            nst=nstim;
        end
        if bayes_opt
            bopars = tapas_fitModel_CAB([], u(1:nst,:), prc_model, obs_model, opt_algo); %BAYES OPTIMAL
        else
            bopars = tapas_fitModel_CAB(y(1:nst,:), u(1:nst,:), prc_model, obs_model, opt_algo);
        end
        bopars.conds=condi;
        save(fullfile(aname,sname),'bopars');
        % PLOTS
        %tapas_hgf_binary_plotTraj(bopars)
        %tapas_fit_plotResidualDiagnostics(bopars);
    end
    
    % PLOTS
    
    %NOTES
    % softmax doesn't like multiple columns in U
    

    %alpha = 10; %uncertain
    %alpha = 0.5; %a bit uncertain
    %alpha = 0.05; %certain
    %sim = tapas_simModel(u(2:end), 'tapas_hgf', [-0.0197 1 1 0.09 1 1 0 0 0 1 1 -1.1537 -6.2314 -6.2823 alpha], 'tapas_softmax_binary',[log(48)]);
    %tapas_hgf_plotTraj(sim)


    % change the responses simulated to make 0/1 changes less perceptible
    %prc_model = 'tapas_hgf_config';
    %obs_model = 'tapas_gaussian_obs_config';
    %obs_model = 'tapas_softmax_binary_config';
    %opt_algo = 'tapas_quasinewton_optim_config';
    %bopars = tapas_fitModel(sim.y, u, prc_model, obs_model, opt_algo);
    %bopars = tapas_fitModel(zeros(length(sim.y),1), u, prc_model, obs_model, opt_algo);
    %tapas_hgf_plotTraj(bopars)
end

if acc
    xlswrite(['accuracy' opt_name '.xlsx'],acc_all);
    xlswrite(['reactiontime' opt_name '.xlsx'],rt_all);
end

if sdt
    xlswrite(['dprime' opt_name '.xlsx'],dp_all);
    xlswrite(['criterion' opt_name '.xlsx'],cri_all);
    xlswrite(['hitrate' opt_name '.xlsx'],hit_all);
    xlswrite(['farate' opt_name '.xlsx'],fa_all);
end

if hgf
    fnames = dir(fullfile(aname,['*_bopars' opt_name '.mat']));
    clear rs % struct 
    clear rc % cell
    for f=1:length(fnames)
        load(fullfile(aname,fnames(f).name));
        rs(f)=bopars;
        rc{f}=bopars;
    end
    tapas_fit_plotCorr_group(rs);
    rc=tapas_bayesian_parameter_average_CAB(1,rc);
    tapas_fit_plotCorr(rc);
end
disp('FINISHED');
% LME
%LMEs can be used to calculate Bayes factors by exponentiating the difference in LME
%between two models applied to the same dataset. For example, an LME difference of 3
%implies a Bayes factor of about 20.
%For a fixed-effects analysis with several datasets (e.g., from different subjects), add up
%the LMEs for the different datasets and compare the LME sums.