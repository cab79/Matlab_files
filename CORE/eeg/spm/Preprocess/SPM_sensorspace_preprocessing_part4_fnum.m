clear all
restoredefaultpath
run('C:\Data\Matlab\Matlab_files\CORE\CORE_addpaths')

%% SPECIFY DATA
filepath = 'C:\Data\CORE\EEG\ana\spm\SPMdata'; 
outpath = 'C:\Data\CORE\EEG\ana\spm\SPMdata\sensorimages'; 
batchpath = 'C:\Data\Matlab\Matlab_files\CORE\EEG\spm\Sensor';

% prefix, middle part, or suffix of files to load (or leave empty) to select a subset of files in
% the folder
fpref = 'spm12_fnums_flip';conds=1:8; 
%fpref = 'spm12_blockmismatch';conds=1:12; 
fmid = '';
fsuff = '4_merged_cleaned.mat';
%fsuff = '_2_merged_cleaned.mat';


%% SPECIFY OPTIONS

% use spm auto artefact rejection
spmart = 0;
% use FT manual artefact rejection
ftart = 0;
% output images?
imout=1;
% transform
arcsinh_data = 0;
arcsinh_images=1;

% output type: 'average' (to create average) 'useaverage' (to use existing averaged file) or 'singletrial'
outputtype = 'useaverage'; % CAREFUL WITH USEAVERAGE: data might not be baselined correctly

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
%timewin = [-200 899]; % empty will include whole epoch
timewin = [-200 299]; % empty will include whole epoch
basewin = [-200 0]; % empty will not baseline correct

%smooth output images (specify FWHM or 0 for no smoothing)
spmsmooth = 12;
delete_unsmoothed = 1;

%% RUN
fname=[fpref '*' fmid  '*' fsuff];
fname=strrep(fname,'**','*');
files = dir(fullfile(filepath,fname));
files_ana=1:length(files);

for f = files_ana
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

for f = files_ana
    fname = files(f).name;
 
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
        
        % prepare file: order conditions
        conditionlabels = cellfun(@num2str, num2cell(conds), 'UniformOutput', false);
        if ~isempty(conditionlabels)
            matlabbatch{1}.spm.meeg.preproc.prepare.D = {S.D};
            matlabbatch{1}.spm.meeg.preproc.prepare.task{1}.sortconditions.label = conditionlabels;
            try
                spm_jobman('initcfg')
                spm_jobman('run',matlabbatch);
            end
        end
    elseif strcmp(outputtype,'useaverage')
        clear S;
        S.prefix = 'm';
        fname = [S.prefix fname];
    end
   
    
    % Baseline Correction
    clear S D
    S.D = fullfile(filepath,fname);
    S.timewin = basewin;
    S.save = 0; % save in separate file
    S.prefix = 'b'; % for output, only if saving a new file
    spm_eeg_bc(S);
    
    % transform data
    if arcsinh_data
        D = spm_eeg_load(fullfile(filepath,fname));
        D(:) = log(D(:)+sqrt(D(:).^2+1));
        save(D)
    end
    
    if imout

        % Convert to 3D images (2 space and 1 time dimension)
        %==========================================================================
        clear S D;
        S.mode = mode;
        S.prefix = '';
        if ~isempty(freqwin); 
            S.freqwin = freqwin; 
            S.prefix = [S.prefix 'f' num2str(freqwin(1)) '_' num2str(freqwin(2)) '_'];
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
        if ~any(strcmp(outputtype,{'useaverage','average'}))
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
                matlabbatch{1,1}.spm.spatial.smooth.data = {fullfile(outpath,imgfold,[imgfiles{nf} ',1'])};
                spm_jobman('initcfg')
                spm_jobman('run',matlabbatch);
                if delete_unsmoothed 
                    delete(fullfile(outpath,imgfold,imgfiles{nf}));
                end
            end
        end 
        
        
        % transform images
        if arcsinh_images
            d = dir(fullfile(outpath,imgfold,'*condition*'));
            imgfiles = {d.name}';

            for nf = 1:length(imgfiles)
                imfile = fullfile(outpath,imgfold,imgfiles{nf});
                spm_imcalc_ui(imfile,imfile,'i1+sqrt(i1.^2+1)');
            end
        end
        
    end
end
