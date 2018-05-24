clear all
dbstop if error
%% SPECIFY DATA
filepath = 'C:\Data\Catastrophising study\SPMdata'; 
outpath = 'C:\Data\Catastrophising study\SPMdata\sensorimages_nosmooth'; 
batchpath = 'C:\Data\Matlab\Matlab_files\Cata study\SPManalysis\Sensor';

% prefix, middle part, or suffix of files to load (or leave empty) to select a subset of files in
% the folder
fpref = 'spm12_csd_';
fmid = '';
fsuff = 'cleaned_SPNall.mat';
%fsuff = 'cleaned_trialNmatch.mat';

%% SPECIFY OPTIONS

% use spm auto artefact rejection
spmart = 0;
% use FT manual artefact rejection
ftart = 0;
% overwrite previous versions?
%overwrite = 0;

% output type: 'average' (to create average) 'useaverage' (to use existing averaged file), or 'singletrial'
outputtype = 'average'; % CAREFUL WITH USEAVERAGE: data might not be baselined correctly

% mode - type of images to generate. One of:
%                'scalp x time'
%                'scalp x frequency' (average over time)
%                'scalp' (average over time and frequency)
%                'source' (average over time and frequency)
%                'time x frequency' (average over channels)
%                'time' (1D average over channels, frequency)
%                'frequency' (1D average over channels, time)
%                'average' (average over all dimensions to get a single
%                           number)
mode = 'scalp x time';

% time and frequecy windows
freqwin = []; % empty if not requiring freq analysis
freqres = 0; % freq resolution
%timewin = [-5500 -2500]; % empty will include whole epoch
%basewin = [-5500 -5000]; % empty will not baseline correct
timewin = [-3000 -2]; % empty will include whole epoch
basewin = [-3000 -2500]; % empty will not baseline correct
%timewin = [-5500 1500]; % empty will include whole epoch
%basewin = [-5500 -5000]; % empty will not baseline correct

%smooth output images (specify FWHM or 0 for no smoothing)
spmsmooth = 1;
delete_unsmoothed = 0;

%% RUN
files = dir(fullfile(filepath,[fpref '*' fmid  '*' fsuff]));

for f = 1:length(files)
    fname = files(f).name;
    if spmart
        
        % Artefact rejection parameters
        %==========================================================================
        clear S;
        S.mode = 'mark'; % or 'reject'
        S.badchanthresh = 0.2;
        S.methods(1).channels = 'all';
        S.methods(1).fun = 'threshchan';
        S.methods(1).settings.threshold = 100;
        S.methods(1).settings.excwin = 1000;
        S.methods(2).channels = 'all';
        S.methods(2).fun = 'jump';
        S.methods(2).settings.threshold = 100;
        S.methods(2).settings.excwin = 1000;
        %S.methods(3).channels = 'all';
        %S.methods(3).fun = 'peak2peak';
        %S.methods(3).settings.threshold = 200;
        S.append = false;
        S.prefix = '';
        
        
        % Artefact rejection 
        %==========================================================================
        S.D = fullfile(filepath,fname);
        fprintf('Artefact rejection: subject %i\n',f);
        D = spm_eeg_artefact(S);
        files(f).name = [S.prefix files(f).name];
    end
    if ftart
        D = spm_eeg_load(fullfile(filepath,fname));
        D = FTrejman(D,[0 0],'SPM');
        save(D)
    end
end

for f = 1:length(files)
    fname = files(f).name;
    
    if ~isempty(freqwin)
        S.D                = fullfile(filepath,fname); %- MEEG object or filename of M/EEG mat-file with
        S.channels         = 'All'; %- cell array of channel names. Can include generic
        try
            S.frequencies      = freqwin(1):freqres:freqwin(2); %- vector of frequencies of interest
        catch
            S.frequencies      = freqwin;
        end
        S.timewin          = timewin; %- time window of interest in PST in ms.
        S.method           ='morlet'; %- name for the spectral estimation to use. This
            %S.settings.subsample   =1;%- factor by which to subsample the time axis (default - 1)
            %S.settings.freqres     =freqres/2;%- frequency resolutions (plus-minus for each frequency
            %S.settings.order       =3;%- butterworth filter order (can be a vector with a value
        S.phase            =0; %- also save phase dataset (1) or not (0)
        prefix           ='tf_'; %- prefix added before the standard prefix (tf_ or tph_)
        S.prefix           =''; %- prefix added before the standard prefix (tf_ or tph_)
        spm_eeg_tf(S)
        fname = [prefix fname];
        
        % Baseline correct
        S.D                =fullfile(filepath,fname);%- MEEG object or filename of M/EEG mat-file
        S.method           ='Log';%- 'LogR', 'Diff', 'Rel', 'Log', 'Sqrt', 'None'
        S.timewin          =basewin;%- 2-element vector: start and stop of baseline (ms)
        S.pooledbaseline   =1;%- take the baseline individually for each trial: doi: 10.1111/ejn.13179
        S.prefix           ='r';%- prefix for the output file (default - 'r')
        spm_eeg_tf_rescale(S)
        fname = [S.prefix fname];
        
    end
 
    if strcmp(outputtype,'average')

        % Average over trials per condition
        %==========================================================================
        clear S D;
        S.prefix = 'm';
        S.review = 0;
        
        % turn these on for robust averaging
        %S.robust.ks = 3;
        %S.robust.bycondition = true;
        %S.robust.savew = false;
        %S.robust.removebad = true;
        %S.prefix = 'mr';
        
        fprintf('Averaging: subject %i\n',f);
        S.D = fullfile(filepath,fname);
        D = spm_eeg_average(S);
        fname = [S.prefix fname];
        if isfield(S,'robust')
            S.D = fullfile(filepath,fname);
            S.band = 'low';
            S.freq = 50;
            S.prefix='';
            spm_eeg_filter(S);
        end
    elseif strcmp(outputtype,'useaverage')
        clear S;
        S.prefix = 'm';
        fname = [S.prefix fname];
    elseif strcmp(outputtype,'singletrial')
        % do nothing here - single trials with output later
    end
   
    
    % Baseline Correction
    if isempty(freqwin)
        clear S D
        S.D = fullfile(filepath,fname);
        S.timewin = basewin;
        S.save = 0; % save in separate file
        S.prefix = 'b'; % for output, only if saving a new file
        spm_eeg_bc(S);
    end

    % Convert to 3D images (2 space and 1 time dimension)
    %==========================================================================
    clear S D;
    S.mode = mode;
    S.prefix = '';
    if ~isempty(freqwin); 
        S.freqwin = freqwin; 
        try
            S.prefix = [S.prefix 'f' num2str(freqwin(1)) '_' num2str(freqwin(2)) '_'];
        catch
            S.prefix = [S.prefix 'f' num2str(freqwin(1)) '_' ];
        end
    end;
    if ~isempty(timewin); 
        S.timewin = timewin; 
        S.prefix = [S.prefix 't' num2str(timewin(1)) '_' num2str(timewin(2)) '_'];
    end;
    if ~isempty(basewin);
        S.prefix = [S.prefix 'b' num2str(basewin(1)) '_' num2str(basewin(2)) '_'];
    end;
    
    fprintf('Converting to images: subject %i\n',f);
    S.D = fullfile(filepath,fname);
    img_files = spm_eeg_convert2images(S);
    
    % move image folder to outpath
    [~,imgfold,~] = fileparts([S.prefix fname]);
    movefile(fullfile(filepath,imgfold),outpath)
    
    % re-orient the images
    d = dir(fullfile(outpath,imgfold,'condition*'));
    imgfiles = {d.name}';
    for nf = 1:length(imgfiles)
        inputname = fullfile(outpath,imgfold,imgfiles{nf});
        nii=load_nii(inputname);
        
        %permute
        nii.img = permute(nii.img,[2 1 3 4]);
        
        %flip LR
        M = diag(nii.hdr.dime.pixdim(2:5));
        M(1:3,4) = -M(1:3,1:3)*(nii.hdr.hist.originator(1:3)-1)';
        M(1,:) = -1*M(1,:);
        nii.hdr.hist.sform_code = 1;
        nii.hdr.hist.srow_x = M(1,:);
        nii.hdr.hist.srow_y = M(2,:);
        nii.hdr.hist.srow_z = M(3,:);
        
        save_nii(nii,inputname);
    end
        
    % split 4D files into 3D
    if strcmp(outputtype,'singletrial')
        d = dir(fullfile(outpath,imgfold,'condition*'));
        imgfiles = {d.name}';
        for nf = 1:length(imgfiles)
            inputname = fullfile(outpath,imgfold,imgfiles(nf));
            spm_file_split(inputname{:}, fullfile(outpath,imgfold));
            delete(inputname{:});
        end
    end

    if spmsmooth 
        % Smooth images
        %==========================================================================
        load(fullfile(batchpath,'batch_smooth'));
        matlabbatch{1,1}.spm.spatial.smooth.fwhm = [spmsmooth spmsmooth 0]; % space space time
        matlabbatch{1,1}.spm.spatial.smooth.dtype = 0;
        matlabbatch{1,1}.spm.spatial.smooth.im = 1; % implicit mask
        d = dir(fullfile(outpath,imgfold,'condition*'));
        imgfiles = {d.name}';
        
        for nf = 1:length(imgfiles)
            imfile = fullfile(outpath,imgfold,imgfiles{nf});
            %if ~exist(imfile) || overwrite
                matlabbatch{1,1}.spm.spatial.smooth.data = {fullfile(outpath,imgfold,[imgfiles{nf} ',1'])};
                spm_jobman('initcfg')
                spm_jobman('run',matlabbatch);
                if delete_unsmoothed 
                    delete(fullfile(outpath,imgfold,imgfiles{nf}));
                end
            %end
        end
    end 
end
