function SW_connectivity(S)

addpath(genpath('C:\Data\Matlab\osl-core'))
addpath('C:\Data\Matlab\MEG-ROI-nets\ROInets')

% PREREQUISITES: 
% - An SPM.mat file containing a contrast on source-space ERP data 
% - Associated cluster data created from "Extract_clusters_source.m"

%% Requires input structure S containing fields as follows. See Cluster_processing script for examples.
%-------------------------------------------------------------
% S.data_path: root directory in which subject-specific folders are located
% S.spmstats_path: directory in which SPM analysis is saved 
% S.spm_dir: specific folder containing the SPM stats for this analysis
% S.contrasts: contrast name cell array - must match that in Matlabbatch (i.e. from design-batch script). Leave empty to proccess ALL contrasts in Matlabbatch
% S.clus_path: cell array of full paths to clusters to be analysed
%generic cluster data name
gdataname = 'cluster_data*.mat';

%% run

if isempty(S.spm_dir)
    S.spm_paths = dir(fullfile(S.spmstats_path,'*spm*'));
    S.spm_paths = {S.spm_paths(:).name};
elseif ~iscell(S.spm_dir)
    S.spm_paths = {S.spm_dir};
else
    S.spm_paths = S.spm_dir;
end

% For each analysis (time window)
%last_spm_path='';
last_tw=0;
for sp = 1:size(S.spm_paths,1)
    S.spm_path = fullfile(S.spmstats_path,S.spm_paths{sp,1});
    
    % run analysis on combined clusters?
    if strcmp(S.gclusname,'comb_clus.nii')
        anapath = fullfile(S.spmstats_path,['Timewin_' num2str(S.spm_paths{sp,2})]);
        S.contrasts = {'Combined'};
        % continue if this timewindow was just analysis
        if last_tw==S.spm_paths{sp,2}
            continue
        end
        last_tw = S.spm_paths{sp,2};
    else
        anapath = S.spm_path;
    end
    
    %S.spm_path = fullfile(S.spmstats_path,S.spm_paths{sp,1});
    %if strcmp(last_spm_path,S.spm_path)
    %    continue
    %end
    %last_spm_path = S.spm_path;

    load(fullfile(S.spm_path,'SPM.mat'));
    S.Fm = SPM.xX.I; % factor matrix
    S.imglist = SPM.xY.P; % Subject-condition image list
    
    S.clus_path={};
    if isempty(S.contrasts)
        alldir=dir(fullfile(anapath,'*_clusters'));
        for c = 1:length(alldir)
            S.clus_path{c} = fullfile(anapath,alldir(c).name);
        end
    else
        for c = 1:length(S.contrasts)
            S.clus_path{c} = fullfile(anapath,[S.contrasts{c} '_clusters']);
        end
    end

    % For each contrast
    for cldir = 1:length(S.clus_path)
        S.cldir=cldir;

        cdata = dir(fullfile(S.clus_path{cldir},gdataname));

        % For each cluster
        for c = 1:length(cdata)
            disp(['*** Cluster = ' num2str(c) ' / ' num2str(length(cdata))])

            [~,cname,~] = fileparts(cdata(c).name);

            % load
            Sc = load(fullfile(S.clus_path{cldir},cdata(c).name));
            
            % update S with data from Sc
            f = fieldnames(Sc.S);
            for i = 1:length(f)
                if ~isfield(S,f{i})
                    S.(f{i}) = Sc.S.(f{i})
                end
            end
            
            % Get cluster names
            fnames = fieldnames(S.wf);
            
            % only analyse first one for now
            wfall=S.wf.(fnames{1}).wf;
            
            if ~S.ana_singlesub
                % concatenate over subjects and conditions
                wfall = {cat(2,wf{:})};
            end

            orth_wf=cell(size(wfall));
            for w = 1:length(wfall)
                
                wf = wfall{w}';
                
                % get wf size to reshape later
                sizewf = size(wf);
                
                % remove means over time
                wf = demean(wf, 1);

                start_tol=1e-10;
                tol=start_tol;
                while isempty(orth_wf{w})
                    disp(['trying tolerance = ' num2str(tol) ' ...'])
                    % get linearly dependent columns and their indices
                    % ucol is the order of the magnitude of the columns
                    % https://uk.mathworks.com/matlabcentral/answers/108835
                    [ld_wf,ld_sub,ucol]=licols(wf,tol);
                    %uni_ucol = unique(ucol);
                    %for u = uni_ucol'
                    %    mu = mean(wf(:,ucol==u));
                    %    find(ucol==u,1,'first')
                    %end

                    % orthogonalise (can also do this on a sample only to get the
                    % weight matrix)
                    %try
                    switch S.leakageCorrectionMethod
                        case 'symmetric'
                            [orth_wf{w}, ~, ~, W,r,asize] = symmetric_orthogonalise_CAB(ld_wf, 1);
                            disp(['wf' num2str(w) ': rank = ' num2str(r) ', size = ' num2str(asize)])
                        case 'closest' 
                            [L, d, rho, W] = closest_orthogonal_matrix(ld_wf);
                    end
                    %catch
                    if (asize-r)/r > 0.05
                        tol=tol*10;
                    else
                        tol=tol+start_tol;
                    end
                end

                % remove means over time
                orth_wf{w} = demean(orth_wf{w}, 1);
                
                orth_wf{w} =  orth_wf{w}';
            end

            %if islogical(S.ana_singlesub)
            %    if ~S.ana_singlesub
            %        % convert back to subjects/conditions
            %        L=reshape(L',size(L,2),[],sizewf(2));
            %        % convert back to cells per subj/cond
            %        L=squeeze(num2cell(L,[1 2]));
            %    end
            %elseif isnumeric(S.ana_singlesub)
            %    L{w}=L{w};
            %end

            % save 
            S.wf.(fnames{1}).orth_wf = orth_wf;
            S.wf.(fnames{1}).ucol = ucol;
            save(fullfile(S.clus_path{cldir},[cname '.mat']),'S');
            
            
            if 0
                Cnii = load_nii(fullfile(cdir,[fnames{1} '.nii']));

                img = Cnii.img;
                ucol = S.wf.(fnames{1}).ucol;

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
                save_nii(fullfile(cdir,RCnii),'reduced.nii')

                NCnii = Cnii;
                NCnii.img = new_clus_img;
                save_nii(fullfile(cdir,NCnii),'new_clus.nii')
            end
            
            %% select data
            time = S.wf.comb_clus.time;
            subs = unique(S.Fm(:,1+S.subrow));

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


                CorrMats = run_correlation_analysis(swf, swf_fil, S.Regularize)

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
                                                         S.SubjectLevel.contrasts);
                                                             
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
        
        % reformat results - correlationMats is a cell array of frequency bands
        S.nFreqBands = length(S.frequencyBands);
        S.nSessions = length(subs);  
        correlationMats = reformat_results(mats, S);

        %% Subject-level analysis to average over sessions in a fixed-effects manner
        % will be same output as first-level if there is one session per subject
        correlationMats = do_subject_level_glm(correlationMats, S);

        %% Group-level analysis
        % Find whole group means
        if strcmpi(S.paradigm, 'rest'),
            correlationMats = do_group_level_statistics(correlationMats, S);
        end%if

        % Perform group-level GLM
        if ~isempty(S.GroupLevel),
            correlationMats = do_group_level_glm_CAB(correlationMats, S);
        end%if
        
        % save 
        S.wf.(fnames{1}).correlationMats = correlationMats;
        save(fullfile(S.clus_path{cldir},[cname '.mat']),'S');

        if 0
            %% hierarchical cluster analysis
            swf_dist = pdist(swf_fil);
            swf_dist = squareform(swf_dist)
            swf_link = linkage(swf_dist);
            dendrogram(swf_link)
            c = cophenet(swf_link,swf_dist)
            I = inconsistent(swf_link)
        end
           
    end
end

end


%--------------------------------------------------------------------------
function FirstLevel = run_first_level_glm(CorrMats, designMat, contrasts)
%RUN_FIRST_LEVEL_GLM

% input checking
[nTrials, nRegressors] = size(designMat);
nContrasts             = length(contrasts);
[~, nModes, checkMe]   = size(CorrMats.envCorrelation_z);
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
[rho, prho, prhoReg] = deal(zeros(nModes, nModes, nContrasts));

% run GLM on each edge
for i = 1:nModes,
    for j = i+1:nModes,
        rho(i,j,:) = glm_fast_for_meg(squeeze(CorrMats.envCorrelation_z(i,j,:)), ...
                                      designMat, invXtX, pinvX, useContrasts, 0);
        prho(i,j,:) = glm_fast_for_meg(squeeze(CorrMats.envPartialCorrelation_z(i,j,:)), ...
                                      designMat, invXtX, pinvX, useContrasts, 0);
								  
	    % fill in uninformative values with NaN.
		if hasBadEVs,
			rho(i,j,badContrasts)  = NaN;
			prho(i,j,badContrasts) = NaN;
		end%if
        if isfield(CorrMats, 'envPartialCorrelationRegularized_z'),
        prhoReg(i,j,1:nContrasts) = glm_fast_for_meg(squeeze(CorrMats.envPartialCorrelationRegularized_z(i,j,:)), ...
                                      designMat, invXtX, pinvX, useContrasts, 0); 
        prhoReg(i,j,badContrasts) = NaN;
        else
            prhoReg(i,j) = 0;
        end%if
    end%for
end%for

% symmetrise and reformat
for iContrast = nContrasts:-1:1,
    FirstLevel(iContrast).cope.correlation                   = rho(:,:,iContrast) + rho(:,:,iContrast)';
    FirstLevel(iContrast).cope.partialCorrelation            = prho(:,:,iContrast) + prho(:,:,iContrast)';
    FirstLevel(iContrast).cope.partialCorrelationRegularized = prhoReg(:,:,iContrast) + prhoReg(:,:,iContrast)';
end%for

end%run_first_level_glm

