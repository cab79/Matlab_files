function DCM = DCM_invert_singlesubject(S)

spm('defaults','EEG');

% Data and analysis directories
%--------------------------------------------------------------------------


try
    Pbase = S.Pbase;
catch
    Pbase = '';
end
try
    Pdata = fullfile(Pbase,S.Pdata);
catch
    Pdata = Pbase;
end
try
    Panalysis = fullfile(Pbase,S.Panalysis);
catch
    Panalysis = Pbase;
end
try
    Pfile = fullfile(Pdata,S.Pfile);
catch
    error('Specify the name of the file as S.Pfile');
end


cd(Panalysis)

% Data filename
%--------------------------------------------------------------------------
DCM.xY.Dfile = Pfile;

% prepare file: order conditions
if isfield(S,'sortconds')
    if ~isempty(S.sortconds)
        matlabbatch{1}.spm.meeg.preproc.prepare.D = {Pfile};
        matlabbatch{1}.spm.meeg.preproc.prepare.task{1}.sortconditions.label = S.sortconds;
        spm_jobman('initcfg')
        spm_jobman('run',matlabbatch);
    end
end

% Parameters and options used for setting up model

%--------------------------------------------------------------------------
%Model type: ’ERP’ is DCM for evoked responses; cross-spectral densities (CSD); induced responses (IND); phase coupling (PHA)
DCM.options.analysis = 'ERP'; % analyze evoked responses
%ERP is standard. CMC (default) and CMM are canonical microcircuit models based on the idea of predictive coding.
DCM.options.model    = 'ERP'; % ERP model
% Select single equivalent current dipole (ECD) for each source, or you use a patch on the 3D cortical surface (IMG).
DCM.options.spatial  = 'ECD'; % spatial model
DCM.options.trials   = [1 2]; % index of ERPs within ERP/ERF file
DCM.options.Tdcm(1)  = 0;     % start of peri-stimulus time to be modelled
DCM.options.Tdcm(2)  = 200;   % end of peri-stimulus time to be modelled
% A projection of the data to a subspace is used to reduce the amount of data, using use the 
% principal components of the prior covariance of the data; select the number of modes you wish to keep. The default is 8.
DCM.options.Nmodes   = 8;     % nr of modes for data selection
%Select 0 for options.h to detrend and just model the mean. 
%Otherwise select the number of discrete cosine transform (DCT) terms you want to use to model low-frequency drifts (> 0). 
%DCT is similar to a fourier transform. 1 = 1 cycle of the cosine wave per time-period; 4 = 1 to 4 cycles.
DCM.options.h        = 1;     % nr of DCT components
% Onset parameter: to avoid attempting to model very early small deflections. 
% changing the onset prior has an effect on how your data are fitted. 
% default value of 60 ms onset time is a good value for many evoked responses where the first large deflection is seen around 100 ms. 
% However, this value is a prior, i.e., the inversion routine can adjust it. 
% Type several numbers (identical or not) only if there were several input stimuli (e.g. visual and auditory) – can be connected to different model sources.
DCM.options.onset    = 60;    % selection of onset (prior mean)
% Duration (sd): makes it possible to vary the width of the input volley, separately for each of the inputs. 
% This can be used to model more closely the actual input structure (e.g. a long stimulus).
DCM.options.dur      = 16;    % Dispersion (sd)
DCM.options.D        = 1;     % downsampling - speeds up computation
DCM.options.han      = 1;     % Hanning removes the effect of beginning and end responses in time window.
% lock ECD orientations by introducing prior correlations. Useful for modelling bilateral symmetric sources (e.g., auditory cortices).
DCM.options.symmetry = 0;
% lock experimental effects by introducing prior correlations.
% ensures that all the changes in connectivity are the same. 
% This is useful when there is a specific hypothesis that some experimental factor increases (or decreases) all connection strengths.
DCM.options.lock = 0;
%“Optimise source locations” only works in combination with the “ECD” option and allows DCM more freedom with moving the dipoles as part of the optimisation process.
DCM.options.location = 0;
% Trial-specific inputs (C parameter estimation)
DCM.options.multiC = 0;

%--------------------------------------------------------------------------
% Data and spatial model
%--------------------------------------------------------------------------
DCM  = spm_dcm_erp_data(DCM);

%--------------------------------------------------------------------------
% Location priors for dipoles
%--------------------------------------------------------------------------
% specify the prior source locations (in mm in MNI coordinates). Can also load the prior locations from a file
DCM.Lpos  = [[-42; -22; 7] [46; -14; 8] [-61; -32; 8] [59; -25; 8] [46; 20; 8]];
% enter the source names (one name in one row). 
DCM.Sname = {'left AI', 'right A1', 'left STG', 'right STG', 'right IFG'};
Nareas    = size(DCM.Lpos,2);

%--------------------------------------------------------------------------
% Spatial model
%--------------------------------------------------------------------------
DCM = spm_dcm_erp_dipfit(DCM);

%--------------------------------------------------------------------------
% Specify connectivity model
%--------------------------------------------------------------------------
cd(Panalysis)

% A matrix is the connections of the first trial/event
% NB: column index corresponds to the source area, and the row index to the
% target area.

% forward connections
DCM.A{1} = zeros(Nareas,Nareas);
DCM.A{1}(3,1) = 1;
DCM.A{1}(4,2) = 1;
DCM.A{1}(5,4) = 1;

% backward connections
DCM.A{2} = zeros(Nareas,Nareas);
DCM.A{2}(1,3) = 1;
DCM.A{2}(2,4) = 1;
DCM.A{2}(4,5) = 1;

% lateral connections
DCM.A{3} = zeros(Nareas,Nareas);
DCM.A{3}(4,3) = 1;
DCM.A{3}(3,4) = 1;

% B matrix: Gain modulations of connection strengths as set in the A-matrices
% Models the difference between the first and the other modelled evoked responses.
% For example, for two evoked responses, DCM explains the first response by
% using the A-matrix only. The 2nd response is modelled by modulating these connections 
% by the weights in the B-matrix.
DCM.B{1} = DCM.A{1} + DCM.A{2};
DCM.B{1}(1,1) = 1;
DCM.B{1}(2,2) = 1;

% C matrix: Inputs - can be to one or many areas
DCM.C = [1; 1; 0; 0; 0];

%--------------------------------------------------------------------------
% Between trial effects
%--------------------------------------------------------------------------
%‘0 1’ to set the first trial type as the baseline for which to compare the second. 
% Or if there is no clear baseline condition, ‘-1 1’ to set the baseline as the average of the two.
DCM.xU.X = [0; 1];
DCM.xU.name = {'rare'};

%--------------------------------------------------------------------------
% Invert
%--------------------------------------------------------------------------
DCM.name = ['DCM_' DCM.options.analysis '-ana_' DCM.options.model '-model_' DCM.options.spatial '-spatial'];

DCM      = spm_dcm_erp(DCM);
