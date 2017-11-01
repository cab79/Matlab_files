function mats = source_correlation_analysis(S,subs,time,orth_wf,select_regions)

mats = {};
for s = 1:length(subs)

    % list index of subjects' image files (same index as S.wf)
    ind = find(S.Fm(:,1+S.subrow)==subs(s) & S.Fm(:,1+S.timerow)==S.timelev);
    % individual subject design matrix
    dmat = S.Fm(ind,:);

    % select subjects data and reshape to 3D
    swf = orth_wf(ind);
    swf = cat(3,swf{:});

    % low-pass filtering
    swf_fil = [];
    if 1
        for tr = 1:size(swf,3)
            [swf_fil(:,:,tr), time_ds] = envelope_data_filteronly(swf(:,:,tr),        ...    
                                                   time,            ...
                                                   S.EnvelopeParams);
        end
    else
        swf_fil=swf;
        time_ds=time;
    end


    CorrMats = run_correlation_analysis(swf, swf_fil, S.Regularize);
    
    % Use an empirical null to enable conversion to z-stats
    transformSurrogates = 0; % don't use if data has negative values, e.g. ERPs
    RegParams           = struct('do', 1, ...
                                 'rho', CorrMats.Regularization.mean);
    sigma = find_permutation_H0_distribution_width(swf_fil, ...
                                                   S.nEmpiricalSamples, ...
                                                   RegParams, ...
                                                   transformSurrogates);
    CorrMats.H0Sigma = sigma;

    % conversion of correlations to z-stats
    fprintf(' Converting correlations to normal z-stats\n');
    CorrMats = convert_correlations_to_normal_variables(CorrMats, ...
                                                       sigma,      ...
                                                       S.Regularize.do);


    % use an OSL function to generate the subject-specific design matrix
    trialID = 1:size(dmat,1);
    designMat = oat_setup_designmatrix(struct('Xsummary', {S.SubjectLevel.designSummary}, ...
                                              'trialtypes', trialID));
    % Run first-level GLM
    % Parameters, e.g. corralation ooefficients, are subjected
    % to the first-level design matrix to produce parameters
    % per contrast (COPE) (e.g. summing/averaging over levels of the design matrix))
    CorrMats.firstLevel = run_first_level_glm(CorrMats,    ...
                                             designMat,          ...
                                             S.SubjectLevel.contrasts, S.SubjectLevel.interaction);

    % visualise correlations
    %A=mean(CorrMats.envPartialCorrelation,3);
    %A(logical(eye(size(A)))) = 0; % make diag zero to visualise small numbers
    %figure;imagesc(A)
    if 0
        A=mean(CorrMats.envPartialCorrelationRegularized_z,3);
        A(logical(eye(size(A)))) = 0; % make diag zero to visualise small numbers
        figure;imagesc(A)
        A=std(CorrMats.envPartialCorrelationRegularized_z,[],3);
        A(logical(eye(size(A)))) = 0; % make diag zero to visualise small numbers
        figure;imagesc(A)
    end
    mats{s}{1} = CorrMats; % second cell is frequency index
    mats{s}{1}.sessionName = num2str(s);
end

end
%--------------------------------------------------------------------------
function FirstLevel = run_first_level_glm(CorrMats, designMat, contrasts,interaction)
%RUN_FIRST_LEVEL_GLM

% input checking
[nTrials, nRegressors] = size(designMat);
nContrasts             = length(contrasts); nInter = size(interaction,1);
[nModes1, nModes2, checkMe]   = size(CorrMats.envCorrelation_z);
assert(checkMe == nTrials,         ...
      [mfilename ':LostTrials'],   ...
      'Number of trials must match columns of design matrix. \n');

assert(iscell(contrasts),               ...
       [mfilename ':NonCellContrasts'], ...
       'Contrasts must be a cell array. \n');
assert(all(cellfun(@length, contrasts) == nRegressors), ...
       [mfilename ':BadContrastFormat'],                ...
       'All contrasts must have the same length as number of regressors. \n');
   
% make sure contrasts are formatted as a cell array of column vectors
useContrasts = cell(1,nContrasts);
for iContrast = 1:nContrasts,
    useContrasts{iContrast} = contrasts{iContrast}(:);
end%for

% create interactions
useInter = cell(1,nInter);
for iInter = 1:nInter,
    highlevel = contrasts{interaction(iInter,1)}(:);
    uhl = unique(highlevel);
    for hl = 1:length(uhl)
        hl_ind = find(highlevel==uhl(hl));
        lowlevel = zeros(length(contrasts{interaction(iInter,2)}),1);
        lowlevel(hl_ind) = contrasts{interaction(iInter,2)}(hl_ind);
        useInter{iInter,hl} = lowlevel;
    end
end%for
   
% Precompute some helpful things
XtX       = designMat' * designMat;
[RXtX, p] = chol(XtX);
if ~p,
	invXtX = RXtX \ (RXtX' \ eye(nRegressors));
	pinvX  = RXtX \ (RXtX' \ designMat');
	hasBadEVs    = false;
	badContrasts = false(nContrasts, 1);
else
	% design matrix was rank deficient
	% is that because we have missing information for certain trial types?
	badEVs    = all(0 == designMat);
	hasBadEVs = any(badEVs);
	if hasBadEVs,
		warning([mfilename ':MissingTrials'],                   ...
			    '%s: file %s is missing trials for %d EVs. \n', ...
				mfilename, fileName, sum(badEVs));
		badContrasts = logical(cellfun(@(C) any(C(badEVs)), useContrasts));
		invXtX = pinv(XtX);
		pinvX  = invXtX * designMat';
	else
		error([mfilename ':RankDeficientDesign'],                     ...
			  ['%s: the design matrix is rank deficient. ',           ...
			   'Check that you''ve specified your EVs sensibly. \n'], ...
			  mfilename);
	end%if
end%if

% declare memory
[rho, prho, prhoReg] = deal(zeros(nModes1, nModes2, nContrasts+length(useInter)));

% run GLM on each edge
for i = 1:nModes1,
    for j = i+1:nModes2,
        rho(i,j,1:nContrasts) = glm_fast_for_meg(squeeze(CorrMats.envCorrelation_z(i,j,:)), ...
                                      designMat, invXtX, pinvX, useContrasts, 0);
        prho(i,j,1:nContrasts) = glm_fast_for_meg(squeeze(CorrMats.envPartialCorrelation_z(i,j,:)), ...
                                      designMat, invXtX, pinvX, useContrasts, 0);
                                  
        rho(i,j,nContrasts+1:nContrasts+length(useInter)) = glm_fast_for_meg(squeeze(CorrMats.envCorrelation_z(i,j,:)), ...
                                      designMat, invXtX, pinvX, useInter, 0);
        prho(i,j,nContrasts+1:nContrasts+length(useInter)) = glm_fast_for_meg(squeeze(CorrMats.envPartialCorrelation_z(i,j,:)), ...
                                      designMat, invXtX, pinvX, useInter, 0);
                                  
								  
	    % fill in uninformative values with NaN.
		if hasBadEVs,
			rho(i,j,badContrasts)  = NaN;
			prho(i,j,badContrasts) = NaN;
		end%if
        if isfield(CorrMats, 'envPartialCorrelationRegularized_z'),
        prhoReg(i,j,1:nContrasts) = glm_fast_for_meg(squeeze(CorrMats.envPartialCorrelationRegularized_z(i,j,:)), ...
                                      designMat, invXtX, pinvX, useContrasts, 0); 
        prhoReg(i,j,nContrasts+1:nContrasts+length(useInter)) = glm_fast_for_meg(squeeze(CorrMats.envPartialCorrelationRegularized_z(i,j,:)), ...
                                      designMat, invXtX, pinvX, useInter, 0); 
        prhoReg(i,j,badContrasts) = NaN;
        else
            prhoReg(i,j) = 0;
        end%if
    end%for
end%for

% symmetrise and reformat
for iContrast = nContrasts+length(useInter):-1:1,
    FirstLevel(iContrast).cope.correlation                   = rho(:,:,iContrast) + rho(:,:,iContrast)';
    FirstLevel(iContrast).cope.partialCorrelation            = prho(:,:,iContrast) + prho(:,:,iContrast)';
    FirstLevel(iContrast).cope.partialCorrelationRegularized = prhoReg(:,:,iContrast) + prhoReg(:,:,iContrast)';
end%for

end%run_first_level_glm

