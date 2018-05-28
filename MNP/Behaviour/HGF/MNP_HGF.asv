function MNP_HGF(D,S)

dbstop if error
%addpath('C:\Data\Matlab\Matlab_files\CORE\Experiment');
addpath('C:\Data\Matlab\HGF\HGFv5.0');
addpath('C:\Data\Matlab\Matlab_files\MNP\Behavioural\HGF\HGF_functions')

dname='C:\Data\MNP\Pilots\NLTv2\processed';
aname='C:\Data\MNP\Pilots\NLTv2\HGF\3lev_test'; bayes_opt=0; 

overwrite=1;
hgf=1;
sdt=0;
acc = 0;

prc_config = 'GBM_config_3lev'; obs_config = 'logrt_softmax_binary_softmu0_config'; nstim=[];

% load .xlsx file containing 'Participant_ID', 'Group', and covariates
pdatfile = 'C:\Data\MNP\Pilots\Participant_Data.xlsx';
% names of headers in the above xls file:
    subhead = 'Subject';
    grphead = 'Group';
    inchead = 'Include';
% which codes to analyse in 'Include' columns in participant data file?
include_codes = [1];

%% RUN

if ~exist(aname,'dir') && hgf
    mkdir(aname)
end

files_ana = 1;
for f = files_ana'
    
    % stimulus - binary (0 and 1)
%     u=D(1).Sequence.signal(2,:); % second stimulus
%     u(u==1)=0;
%     u(u==2)=1;
%     u(2,:)=1;
%     u=u';
%     y = D(1).Processed.presssignal; % BINARY response
%     y(y==1)=1;
%     y(y==2)=0;
%     y=y';

%     % associative learning: u indicates pairings
%     sig=D(1).Sequence.signal(1:2,:); 
%     u(sig(1,:)==sig(2,:)) = 0;
%     u(sig(1,:)~=sig(2,:)) = 1;
%     u(2,:)=1;
%     u=u';
%     % associative learning - responses indicate pairings
%     y=[]
%     ysig=D(1).Processed.presssignal;
%     y(ysig(1,:)==sig(1,:)) = 1;
%     y(ysig(1,:)~=sig(1,:)) = 0;
%     y(isnan(ysig))=nan;
%     y=y';
    
    % associative learning: u indicates outcome and cues
    sig=D(1).Sequence.signal(1:2,:); 
    u(1,sig(2,:)==1) = 0; % outcome
    u(1,sig(2,:)==2) = 1; % outcome
    u(2,:) = sig(2,:); % outcome types
    u(3,:) = sig(1,:); % cues
    u=u';
    y = D(1).Processed.presssignal; % BINARY response
    y(y==1)=1;
    y(y==2)=0;
    y=y';

    %% HGF
    % prc: perceptual; obs:observation; opt:optimisation
    prc_model = prc_config;
    %prc_model = 'tapas_hgf_binary_pu_config';
    obs_model = obs_config;
    %obs_model = 'tapas_bayes_optimal_binary_config'; %BAYES OPTIMAL
    opt_algo = 'tapas_quasinewton_optim_config';
    
    sname = datestr(now,30);
    sname = [S(1).select.subjects{1, 1} '_' sname '_bopars.mat'];
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
        %bopars.conds=condi;
        save(fullfile(aname,sname),'bopars');
        % PLOTS
        tapas_hgf_binary_plotTraj(bopars)
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

% if acc
%     xlswrite(['accuracy' opt_name '.xlsx'],acc_all);
%     xlswrite(['reactiontime' opt_name '.xlsx'],rt_all);
% end
% 
% if sdt
%     xlswrite(['dprime' opt_name '.xlsx'],dp_all);
%     xlswrite(['criterion' opt_name '.xlsx'],cri_all);
%     xlswrite(['hitrate' opt_name '.xlsx'],hit_all);
%     xlswrite(['farate' opt_name '.xlsx'],fa_all);
% end
% 
% if hgf
%     fnames = dir(fullfile(aname,['*_bopars' opt_name '.mat']));
%     clear rs % struct 
%     clear rc % cell
%     for f=1:length(fnames)
%         load(fullfile(aname,fnames(f).name));
%         rs(f)=bopars;
%         rc{f}=bopars;
%     end
%     tapas_fit_plotCorr_group(rs);
%     rc=tapas_bayesian_parameter_average_CAB(1,rc);
%     tapas_fit_plotCorr(rc);
% end
disp('FINISHED');
% LME
%LMEs can be used to calculate Bayes factors by exponentiating the difference in LME
%between two models applied to the same dataset. For example, an LME difference of 3
%implies a Bayes factor of about 20.
%For a fixed-effects analysis with several datasets (e.g., from different subjects), add up
%the LMEs for the different datasets and compare the LME sums.