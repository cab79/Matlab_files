% Fit a DCM to a spectra

clear all
config.datafile        = 'L:\myfolder\mydatafile'; %%% THIS WOULD CORRESPOND TO SOURCE EXTRACTED data file
config.Outputfile      = 'L:\myfolder\DCM_outputfilename'; %%% CAN NAME THis anything

load(config.datafile)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
DCM.name    = config.Outputfile; %% or whatever you like

%% Define the relevant parameters for the DCM
DCM.xY              = [];
DCM.xY.Dfile        = config.datafile;  % data file
DCM.xY.modality     = 'LFP';
DCM.xY.Ic           = [1:8];  %% THERE were 8 channels in this file - you may want to extract 2 e.g. PCC and ACC
DCM.Sname           = {'1' '2' '3' '4' '5' '6' '7' '8'};       %% These are strings that anme the data series
DCM.options.Nmodes  = 8;
DCM.M.dipfit.Ns     = 8;
DCM.M.dipfit.Nc     = 8;

DCM.options.trials   = [1 2 3];        % trial codes - I had three conditions here
DCM.options.analysis = 'CSD';      % CSD stands for cross spectral densities - to fit the average spectrum in the time seris
DCM.options.model    = 'CMC';       % CAN BE CHANGED  LFP or ERP, these latter are quite stable.
DCM.options.spatial  = 'LFP';       % this option is for when you are 
DCM.options.Fdcm     = [2 60];        % Frequency window
DCM.options.Tdcm     = [1 1000];   % Peri-stimulus time window


DCM.options.h        = 1;
DCM.options.D        = 8;        % Down-sampling - this will be done automatically in the DCM files
DCM.options.han      = 0;        % No hanning
DCM.options.loc      = 0;
DCM.options.location = 0;
DCM.options.symmetry = 0;
DCM.xU               = [];
DCM.xU.X             = [0 1 0; 0 0 1]';  %% Design matrix for a three level design 


%% Connectivity - forward up the sequence backward back down the sequence
DCM.A{1}    = [0 0 0 0 0 0 0 0;
               0 0 0 0 0 0 0 0;
               0 0 0 0 0 0 0 0;
               0 0 0 0 0 0 0 0;
               0 0 0 0 0 0 0 0;
               0 0 0 0 0 0 0 0;
               1 1 1 1 1 1 0 0;
               1 1 1 1 1 1 0 0];         %% Forward
DCM.A{2}    = DCM.A{1}'        ;         %% Backward
DCM.A{3}    = zeros(8,8) ;               %% Lateral
DCM.B{1}    = DCM.A{1} + DCM.A{2};            %%Modulatory Effects (experimental differences)
DCM.B{2}    = DCM.A{1} + DCM.A{2};       
DCM.C       = zeros(1,8)';                               %% No explicit inputs in CSD 

%% Invert!!
DCM = spm_dcm_csd(DCM);

%% Save
save(DCM.name,'DCM');
