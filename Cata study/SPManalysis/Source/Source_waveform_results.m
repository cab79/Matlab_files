%% take output from othogonisation and produce reduced nii image
load('C:\Data\Catastrophising study\SPMstats\Source\All_timewin\Combined_clusters\cluster_data1.mat')
clusname = 'comb_clus';
Cnii = load_nii(['C:\Data\Catastrophising study\SPMstats\Source\All_timewin\Combined_clusters\' clusname '.nii']);

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

%% hierarchical cluster analysis

% concatenate over subjects and conditions
%wf = cat(2,S.wf.comb_clus.wf{:});
if S.ana_singlesub
    wf = S.wf.comb_clus.orth_wf;
else
    wf = cat(2,S.wf.comb_clus.orth_wf{:});
end
wf_dist = pdist(wf);
wf_dist = squareform(wf_dist)
wf_link = linkage(wf_dist);
dendrogram(wf_link)
c = cophenet(wf_link,wf_dist)
I = inconsistent(wf_link)

%% correlation analysis
addpath('C:\Data\Matlab\MEG-ROI-nets\ROInets')
Regularize.do = 1;
Regularize.method = 'Bayesian';
Regularize.Prior.a = 1/3;
Regularize.Prior.b = 0;
CorrMats = run_correlation_analysis(wf, wf, Regularize)
A=CorrMats.envPartialCorrelation;
A(logical(eye(size(A)))) = 0; % make diag zero to visualise small numbers
figure;imagesc(A)

% Use an empirical null to enable conversion to z-stats
transformSurrogates = 0; % don't use if data has negative values, e.g. ERPs
RegParams           = struct('do', 1, ...
                             'rho', CorrMats.Regularization.mean);
nEmpiricalSamples = 8;
sigma = find_permutation_H0_distribution_width(wf, ...
                                               nEmpiricalSamples, ...
                                               RegParams, ...
                                               transformSurrogates);
CorrMats.H0Sigma = sigma;
    
%% conversion of correlations to z-stats
fprintf(' Converting correlations to normal z-stats\n');
CorrMats = convert_correlations_to_normal_variables(CorrMats, ...
                                                   sigma,      ...
                                                   Regularize.do);
A=CorrMats.envPartialCorrelation_z;
A(logical(eye(size(A)))) = 0; % make diag zero to visualise small numbers
figure;imagesc(A)
                                                               