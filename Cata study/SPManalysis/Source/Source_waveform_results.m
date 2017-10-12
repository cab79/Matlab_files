%% take output from othogonisation and produce reduced nii image
close all
clear all
addpath(genpath('C:\Data\Matlab\osl-core'))
cdir = 'C:\Data\Catastrophising study\SPMstats\Source\Timewin_3\Combined_clusters';
cd(cdir)
load(fullfile(cdir,'cluster_data1.mat'));

if 0
    clusname = 'comb_clus';
    Cnii = load_nii(fullfile(cdir,[clusname '.nii']));

    img = Cnii.img;
    ucol = S.wf.(clusname).ucol;

    [uimg,iU,iI] = unique(img);
    reduced_img = zeros(size(img));
    new_clus_img = zeros(size(img));
    for u = 1:length(uimg)
        if uimg(u)==0
            continue
        else
            reduced_img(iI==u) = ucol(uimg(u));
            % when does ucol have multiple repetitions of it's value?
            if sum(ucol==ucol(uimg(u)))>1
                new_clus_img(iI==u) = ucol(uimg(u));
            end
        end
    end
    RCnii = Cnii;
    RCnii.img = reduced_img;
    save_nii(RCnii,'reduced.nii')

    NCnii = Cnii;
    NCnii.img = new_clus_img;
    save_nii(NCnii,'new_clus.nii')
end

%% select data
wf = S.wf.comb_clus.orth_wf;
time = S.wf.comb_clus.time;
subs = unique(S.Fm(:,1+S.subrow));


    
for w = 1:length(wf)

    
    %% low-pass filtering
    if 1
    Settings.EnvelopeParams.useEnv = false; % envelope the data?
    Settings.EnvelopeParams.useFilter = true; % use a more sophisticated filter than a sliding window average
    Settings.EnvelopeParams.takeLogs  = false; % results in complex numbers with ERPs!
    Settings.EnvelopeParams.windowLength  = 0.2; % in seconds; 1/wl = maxfreq
    [wf_fil, time_ds] = envelope_data_filteronly(wf,        ...    
                                           time,            ...
                                           Settings.EnvelopeParams);
    else
        wf_fil=wf;
        time_ds=time;
    end


    %% correlation analysis
    addpath('C:\Data\Matlab\MEG-ROI-nets\ROInets')
    Regularize.do = 1;
    %Regularize.method = 'Bayesian';
    Regularize.method = 'Friedman';                     % Regularization approach to take. {'Friedman' or 'Bayesian'}
    Regularize.Prior.a = 1/3;
    Regularize.Prior.b = 0;
    Regularize.path          = 0.001;                          % This specifies a single, or vector, of possible rho-parameters controlling the strength of regularization. 
    Regularize.adaptivePath  = false;                          % adapt the regularization path if necessary
    CorrMats = run_correlation_analysis(wf, wf_fil, Regularize)
    A=CorrMats.envPartialCorrelation;
    A(logical(eye(size(A)))) = 0; % make diag zero to visualise small numbers
    figure;imagesc(A)

    % Use an empirical null to enable conversion to z-stats
    transformSurrogates = 0; % don't use if data has negative values, e.g. ERPs
    RegParams           = struct('do', 1, ...
                                 'rho', CorrMats.Regularization.mean);
    nEmpiricalSamples = 8;
    sigma = find_permutation_H0_distribution_width(wf_fil, ...
                                                   nEmpiricalSamples, ...
                                                   RegParams, ...
                                                   transformSurrogates);
    CorrMats.H0Sigma = sigma;

    %% conversion of correlations to z-stats
    fprintf(' Converting correlations to normal z-stats\n');
    CorrMats = convert_correlations_to_normal_variables(CorrMats, ...
                                                       sigma,      ...
                                                       Regularize.do);
    A=CorrMats.envPartialCorrelationRegularized_z;
    A(logical(eye(size(A)))) = 0; % make diag zero to visualise small numbers
    figure;imagesc(A)
end
                         

if 0
    %% hierarchical cluster analysis
    wf_dist = pdist(wf_fil);
    wf_dist = squareform(wf_dist)
    wf_link = linkage(wf_dist);
    dendrogram(wf_link)
    c = cophenet(wf_link,wf_dist)
    I = inconsistent(wf_link)
end