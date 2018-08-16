function [S,D] = CORE_decEEG_import(S,D,ycol)
% imports reconstructed variables from EEG decoding results file; creates
% HGF y variable from column ycol onwards.

file=load(fullfile(S.path.stats,S.file));
stats = file.stats;

for d = 1:length(D)
  
    % get data
    if S.use_group_recons
        recons = stats.biem.alldata(end).grprecons{d,1}; % output of EEG decoding in trial order of conData(testidx)
    else
        recons = stats.biem.alldata(end).recons{d,1}; % output of EEG decoding in trial order of conData(testidx)
    end
    orig = stats.biem.alldata(end).orig{d,1}; % original predictor used to train the EEG encoding model (extracted here for validation purposes)
    tnums = stats.trialinfo{1}.tnums{d,1}; % trial indices remaining in EEG data after cleaning etc., according to original experimental trial order
    idx = stats.trialinfo{1}.idx{d,1}; % index of tnums for both training and testing data
    testidx = stats.trialinfo{1}.testidx{d,1}; % index of idx containing the testing trials
    
    % re-order
    [~,ord] = sort(idx(testidx));
    recons_ord = recons(ord);
    orig_ord = orig(ord);

    % construct y (HGF response variable)
    ycols = ycol:ycol-1+size(recons_ord,2);
    nk = length(D(d).HGF.u);
    if ~isfield(D(d).HGF,'y')
        D(d).HGF.y = nan(nk,ycols);
    else
        D(d).HGF.y(1:nk,ycols) = nan;
    end
    D(d).HGF.y(tnums,ycols)=recons_ord;
    
    % validation: only use if recon is mismatch
    if 0
        % for validation purposes only
        orig_y = nan(nk,1);
        orig_y(tnums,1) = orig_ord;
        curr_u = D(d).HGF.u(:,1);
        curr_u_nonan = curr_u; curr_u_nonan(isnan(orig_y))=[];
        orig_y_nonan = orig_y; orig_y_nonan(isnan(orig_y))=[];
        D(d).orig_y=orig_y;

        % plot
        figure(1); clf
        subplot(2,2,1); plot(recons,'b'); hold on; plot(orig,'r'); title('original (red) vs. recons (blue)')
        subplot(2,2,2); plot(D(d).HGF.y); title('recons ordered as per HGF')
        subplot(2,2,3); scatter(D(d).HGF.u(:,1),D(d).HGF.y(:,ycols)); title('should be a mean difference (y axis) between 0 and 1 (x axis)')
        subplot(2,2,4); plot(curr_u_nonan,'b'); hold on; plot(orig_y_nonan,'r'); title('overlap of u with recon - should be no blue showing')
        if any(curr_u_nonan~=round(orig_y_nonan))
            dbstop if error
            find(curr_u_nonan~=orig_y_nonan)
            error('original predictor not consistent with HGF input')
        end
        pause(0.1)
    end
end
try
close(1)
end