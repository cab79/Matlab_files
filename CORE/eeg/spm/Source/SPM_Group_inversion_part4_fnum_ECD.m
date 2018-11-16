%% SPM source analysis: setup
clear all;close all
restoredefaultpath
run('C:\Data\Matlab\Matlab_files\CORE\CORE_addpaths')
dbstop if error
%load('C:\Data\Catastrophising study\meegfid\fidmat2.mat');
%for fm = 1:length(fidmat)

%% generic directories for all analyses for this study
%-------------------------------------------------------------
% root directory in which SPM data files are located
S.filepath = 'C:\Data\CORE\EEG\ana\spm\SPMdata'; 
% place to save source images/data
S.outpath = 'C:\Data\CORE\EEG\ana\spm\SPMdata\ECD_outputs'; 
% load .xlsx file containing 'Participant_ID', 'Group', and covariates
S.pdatfile = 'C:\Data\CORE\Participants\Participant_data.xlsx';
%fiducials directory
S.fid_dir='C:\Data\CORE\eeg\ana\spm\SPMdata\meegfid';

%% specific directory and file information for this analysis
%-------------------------------------------------------------
% prefix, middle part, or suffix of files to load (or leave empty) to select a subset of files in
% the folder
S.fpref = 'mspm12_fnums_flip';
S.fmid = '';
%S.fsuff = '4_cleaned_tm.mat';
S.fsuff = {
    '_4_merged_cleaned.mat';
    };

%% specific settings for this analysis
% which codes to analyse in 'Include' columns in participant data file?
S.include_codes = [1];
S.grps = {1;2}; %inversion on each group separately: separate with colon. Otherwise separate with comma
% time and frequecy windows
S.freqwin = [0 100]; % high pass and low pass filter (not applied for ECD)
S.timewin_file = {'C:\Data\CORE\eeg\ana\ERP\lat_results_grandavg_20181110T114702.xlsx','event1_peak1',5};%[30 50]; % empty will include whole epoch. Best to narrow this as much as possible to the range of interest. ECD averages over this window so should be short. Can be a single number.
S.cond = {{'1'},{'2'},{'3'},{'4'},{'5'},{'6'},{'7'},{'8'}}; % ECD: finds separate dipoles for each
%S.cond = {{'8'}};
S.basewin = [-200 0]; % empty will not baseline correct and will not produce a baseline image.
S.baseest = 0; % 1 = baseline estimated in separate model
% latency ranges of images to procedure. Will only produce images within
% the S.timewin range.
S.imout = []; % select which to produce on this run, or leave empty to run all
S.images_out = { % empty for ECD
    {
     }
};

sensorpath = '';

%end

S.sourceprior = 'ECD'; % Priors on sources, e.g. ECD, MSP, GS, COH or IID

% IMAGING (MSP/COH/IID) options
S.smooth_prior = 0.2;
S.smooth = 6;%smooth output images (specify FWHM or 0 for no smoothing)
S.type = 'evoked'; % 'evoked', 'induced' or 'trials' (single trial images)
S.group_inversion = 1;
S.Npriors = 4096; % Number of sparse priors (x 1/2 brain)
S.Han = 0; % apply Hanning window
S.run_forward = 0; Msize  = 2;
S.prior_image = 'C:\Data\CORE\eeg\ana\spm\SPMdata\masks\AAL_postcentral_right_6002.nii';
S.arcsinh_data = 0;
S.log_images = 1;

% ECD options
S.prior(1).prior_loc = [39.4 -27.6 63.3]; % dipole prior location in MNI. Leave empty to have a non-informative prior.
S.prior(1).prior_var = 3*[10.4 8.0 6.1]; % dipole prior location variance in MNI (mm).
S.prior(1).prior_mom = []; % prior moment
S.prior(1).prior_momvar = []; % prior moment variance
S.prior(1).sym = 1; % 1=single dipole, 2 = symmetric pair
S.prior(2).prior_loc = [54 -20 20]; % dipole prior location in MNI. Leave empty to have a non-informative prior.
S.prior(2).prior_var = 3*[10 10 10]; % dipole prior location variance in MNI (mm).
S.prior(2).prior_mom = []; % prior moment
S.prior(2).prior_momvar = []; % prior moment variance
S.prior(2).sym = 2; % 1=single dipole, 2 = symmetric pair
S.Niter = 10;

% run options
tw_run = 1%:size(S.timewin,1);
S.figures_on = 0;

%% RUN
if S.figures_on
    set(0,'DefaultFigureVisible','on')
else
    set(0,'DefaultFigureVisible','off')
end

%S.outpath = [S.outpath num2str(fm)];
if ~exist(S.outpath,'dir')
    mkdir(S.outpath);
end

cd(S.outpath)

[~,~,pdata] = xlsread(S.pdatfile);
grp_col = find(strcmp(pdata(1,:),'Group'));
sub_col = find(strcmp(pdata(1,:),'Subject'));
inc_col = find(strcmp(pdata(1,:),'Include'));

%Change char Grp inputs to numbers
grpdat = pdata(2:end,grp_col);
if isnumeric(grpdat{2,1})
    grptype = unique([grpdat{:}]);
    grptype(grptype==0)=[];
    Ngrp = length(grptype);
else
    grptype = unique(grpdat);
    grptype(isempty(grptype))=[];
    Ngrp = length(grptype);
    for g = 1:Ngrp
        grp_idx = cellfun(@(x) any(strcmp(grptype(g),x)), grpdat, 'UniformOutput', 0);
        grpdat(cell2mat(grp_idx)) = {[g]};
    end
end

% find index of subjects to include in analysis
SubInd = cell(Ngrp,1);
Subs = [];
inc_idx = cellfun(@(x) ismember(x,S.include_codes), pdata(2:end,inc_col), 'UniformOutput', 0);
inc_idx = find(cell2mat(inc_idx));

% find subject indices for each specific group
for g = 1:Ngrp
    grp_idx = find(cellfun(@(x) x==g, grpdat, 'UniformOutput', 1));
    SubInd{g,1} = intersect(inc_idx,grp_idx);
    Nsub(g,1) = length(SubInd{g,1});
end

% keep groups split or combine according to S.grps
SubInd_ana = {};
for g = 1:size(S.grps,1)
    SubInd_ana{g} = cell2mat(SubInd([S.grps{g,:}]));
end
Ngrp_ana = length(SubInd_ana);
Nsub_ana=[];

% cycle though S.timewin
for tw = tw_run
    
    if isfield(S,'timewin_file') && ~isempty(S.timewin_file) % latency specified from a file
        [NUM,TXT,RAW]=xlsread(S.timewin_file{1});
        col = find(strcmp(TXT(1,2:end),S.timewin_file{2}));
        timewin = [NUM(:,col)-S.timewin_file{3}, NUM(:,col)+S.timewin_file{3}];
        timewin_fname = TXT(2:end,1);
        timewin_subs = cellfun(@(x) strsplit(x,'_'),timewin_fname,'UniformOutput',false);
        timewin_subs = vertcat(timewin_subs{:, 1});  
        timewin_subs = timewin_subs(:,1);
    else
        timewin = S.timewin(tw,:);
    end
    
    
    try
        basewin = S.basewin(tw,:);
    catch
        basewin = S.basewin(1,:);
    end
    basecon=0;
    if ~isempty(S.images_out)
        images_out = S.images_out{tw};
    else
        images_out={};
    end
    
    % cycle through groups and subjects
    for g = 1:Ngrp_ana
        subID={};
        files={};
        Nsub_ana(g) = length(SubInd_ana{g});
        for s = 1:Nsub_ana(g)
            subID{s} = pdata{SubInd_ana{g}(s)+1,sub_col};
            if isnumeric(subID{s}); subID{s} = num2str(subID{s}); end
            fname=[S.fpref '*' subID{s} '*' S.fmid  '*' S.fsuff{tw}];
            fname=strrep(fname,'**','*');
            fname = dir(fullfile(S.filepath,fname));
            files{s} = fname.name;
            % Baseline Correction
            clear B D
            B.D = fullfile(S.filepath,files{s});
            B.timewin = basewin;
            B.save = 0; % save in separate file
            B.prefix = 'b'; % for output, only if saving a new file
            spm_eeg_bc(B);
        end

        clear functions D;
        % Load data and set method
        %==========================================================================
        %if g==2; st = 11; else st =1;end
        st=1;
        for i = st:length(files)

            fprintf('checking for previous inversions: subject %i\n',i);

            D{i}     = spm_eeg_load(fullfile(S.filepath,files{i}));

            % clear previous inversions
            %----------------------------------------------------------------------
%             try, D{i}.inv{D{i}.val} = rmfield(D{i}.inv{D{i}.val},'inverse' ); end
%             try, D{i}.inv{D{i}.val} = rmfield(D{i}.inv{D{i}.val},'contrast'); end
            
            try, inv = D{i}.inv(1);end
            try, D{i} = rmfield(D{i},'inv'); end
            try, inv{1} = rmfield(inv{1},'inverse' );end
            try, inv{1} = rmfield(inv{1},'inverse' );end
            D{i}.inv = inv;
            
            D{i}.val             = 1;
            if ~strcmp(S.sourceprior,'ECD')
                D{i}.inv{D{i}.val}.method   = 'Imaging';
            end

            % clear redundant models
            %----------------------------------------------------------------------
            %if D{i}.val==1
            %    D{i}.inv = D{i}.inv(D{i}.val);
            %end

            % save forward model parameters
            %----------------------------------------------------------------------
            save(D{i});

        end

        % Check for existing forward models and consistent Gain matrices
        %--------------------------------------------------------------------------
        Nd = zeros(1,length(st:length(files)));
        for i = st:length(files)
            fprintf('checking for foward models: subject %i\n',i);
            if S.run_forward
                Nd(i) = 0;
            else
                try
                    [L, D{i}] = spm_eeg_lgainmat(D{i});
                    Nd(i) = size(L,2);               % number of dipoles
                catch
                    Nd(i) = 0;
                end
            end
            clear L;
        end

        % use template head model where necessary
        %==========================================================================
        if max(Nd > 1024)
            NS = find(Nd ~= max(Nd));            % subjects requiring forward model
        else
            NS = st:length(files);
        end


        %    for i = NS

              
       %     end

            

        %    for i = NS(2:end)
        %        [mod, list] = modality(D{i}, 1, 1);
        %        if ~ismember(inverse.modality, list)
        %            error([inverse.modality ' modality is missing from ' D{i}.fname]);
        %        end
        %    end

            % and save them (assume trials = types)
            %--------------------------------------------------------------------------
            %for i = NS
            %end

        
        if ~isempty(NS)
            for i = NS(1)
                
                cd(D{i}.path);

                

                % use a template head model and associated meshes
                %======================================================================
                D{i}.inv{D{i}.val}.comment =''; % otherwise D.inv gets wiped
                D{i}.inv{D{i}.val}.date =''; % otherwise D.inv gets wiped
                D{i} = spm_eeg_inv_mesh_ui(D{i}, D{i}.val, 1, Msize);
                
                % Select modality
                %==========================================================================
                inverse.modality  = spm_eeg_modality_ui(D{i}, 1, 1);
                
                D{i}.inv{D{i}.val}.inverse = inverse;

                %if ~exist('D','var'); load tempD; end;

                fprintf('Registering and computing forward model (subject: %i)\n',i);

                % Forward model
                %----------------------------------------------------------------------
                %D{i} = spm_eeg_inv_datareg_ui(D{i}, 1);
                load(fullfile(S.fid_dir,'meegfid.mat'));
                %load(fullfile(S.fid_dir,'meegfid.mat'));
                %fm=2
                %meegfid.pnt(1,3) = meegfid.pnt(1,3)+fidmat(fm,1);
                %meegfid.fid.pnt(1,3) = meegfid.fid.pnt(1,3)+fidmat(fm,1);
                %meegfid.fid.fid.pnt(1,3) = meegfid.fid.fid.pnt(1,3)+fidmat(fm,1);
                %meegfid.pnt(2:3,3) = meegfid.pnt(2:3,3)+fidmat(fm,2);
                %meegfid.fid.pnt(2:3,3) = meegfid.fid.pnt(2:3,3)+fidmat(fm,2);
                %meegfid.fid.fid.pnt(2:3,3) = meegfid.fid.fid.pnt(2:3,3)+fidmat(fm,2);

                load(fullfile(S.fid_dir,'newmrifid.mat'));
                %D{i}.fiducials.pnt(2:3,1) = -D{i}.fiducials.pnt(2:3,1);
                useheadshape = 1;
                D{i} = spm_eeg_inv_datareg_ui(D{i},D{i}.val, meegfid,newmrifid, useheadshape);
                %D{i} = spm_eeg_inv_datareg_ui(D{i},1);
                D{i} = spm_eeg_inv_forward_ui_CAB(D{i}, D{i}.val);
                clear functions
                [L, D{i}] = spm_eeg_lgainmat(D{i});
                clear functions
                % save forward model
                %----------------------------------------------------------------------
                save(D{i});
                clear L;

            end

            for i = NS(2:end)
 
                D{i}.inv{D{i}.val} = D{NS(1)}.inv{D{i}.val};
                %D{i}.inv{1} = rmfield(D{i}.inv{1},'gainmat');

                %if ~exist('D','var'); load tempD; end;

                fprintf('Registering and computing forward model (subject: %i)\n',i);

                %clear functions
                %[L, D{i}] = spm_eeg_lgainmat(D{i});
                %clear functions
                % save forward model
                %----------------------------------------------------------------------
                save(D{i});
                %clear L;

            end
        end
        
        if strcmp(S.sourceprior,'ECD')
            for i = st:length(files)
                if ~isfield(S,'timewin')
                    % index of subject within timewin file
                    subind=[];
                    for fi = 1:length(timewin_subs)
                        if ~isempty(strfind(files{i},timewin_subs{fi}))
                            subind=fi;
                        end
                    end
                    S.timewin = timewin(subind,:);
                    disp(['time window: ' num2str(S.timewin)]);
                end
                D{i} = spm_eeg_inv_vbecd_cab(D{i},S);
            end
            
        else
            %if ~exist('contrasts','var')
                if S.baseest
                    contrasts = {timewin,[]; basewin,[]};
                    for i = st:length(files)
                        D{i}.inv{size(S.timewin,1)+1} = D{i}.inv{tw}; 
                    end
                else
                    contrasts = {timewin,[]};
                end
            %end

            %if ~exist('run_con','var')
                if S.baseest
                    run_con = 1:2;
                else
                    run_con = 1;
                end
            %end

            for con = run_con

                for i = st:length(files)

                    % creates new D.inv for each baseline contrast
                    if con==2
                        D{i}.val=1+size(S.timewin,1);
                    else
                        D{i}.val=1;
                    end

                    woi=[];
                    %D{i}.inv{D{i}.val}.inverse.trials = D{i}.condlist; % Trials
                    D{i}.inv{D{i}.val}.inverse.type   = S.sourceprior;
                    %D{i}.inv{D{i}.val}.inverse.smooth = 0;        % Smoothness of source priors (mm)
                    D{i}.inv{D{i}.val}.inverse.Np     = S.Npriors;
                    D{i}.inv{D{i}.val}.inverse.Han = S.Han; 
                    if ~isempty(S.freqwin)
                        D{i}.inv{D{i}.val}.inverse.hpf = S.freqwin(2); 
                        D{i}.inv{D{i}.val}.inverse.lpf = S.freqwin(1); 
                    end

                    contrast = contrasts{con,1};
                    if  size(contrast,2)==1   
                        if ~isempty(contrasts{con,2})
                            conwin = contrasts{con,2}/2;
                            if size(contrast,1)==1
                                woi              = [contrast-conwin contrast+conwin];
                            else 
                                woi              = [contrast(grpind(i),1)-conwin contrast(grpind(i),1)+conwin];
                            end
                        else 
                            woi              = contrast;
                        end
                    elseif size(contrast,2)==2 && size(contrast,1)>1
                        woi              = [contrast(grpind(i),1) contrast(grpind(i),2)];
                    else
                        woi              = contrast;
                    end
                    woi              = sort(woi);
                    woi     = round([woi(1) woi(end)]);

                    D{i}.inv{D{i}.val}.inverse.woi =woi; %1000*[min(D.time) max(D.time)];

                    D{i}.inv{D{i}.val}.inverse.smooth = S.smooth_prior;

                    % spatial priors
                    if ~isempty(S.prior_image) && i==1
                        Ss.D = D{i};      %- MEEG object or filename of M/EEG mat-file
                        Ss.fmri = S.prior_image;  %- filename of prior (SPM) image to be used
                        Ss.space = 1;  %- native (0) or MNI (1) space (must be same for SPM and GM images)
                        Ss.hthr = 0.5;  %- height threshold of prior image [defaults: 0.5]
                        Ss.ethr = 1;   %- extent threshold of clusters in prior image [default: 1]
                        Ss.ncomp = Inf;  %- maximal number of priors component to be extracted [default: Inf]
                        Ss.disp = 0;  %- whether to display priors on mesh [default: 0]
                        D{i} = spm_eeg_inv_fmripriors(Ss);
                        load(D{i}.inv{D{i}.val}.inverse.fmri.priors);
                        D{i}.inv{D{i}.val}.inverse.pQ = pQ;
                    elseif ~isempty(S.prior_image) 
                        D{i}.inv{D{i}.val}.inverse.pQ = D{1}.inv{D{i}.val}.inverse.pQ;
                    end
                end

            end

            % specify conditions from S.cond
            len_inv = length(D{1}.inv);
            val = 0;
            for cn = 1:length(S.cond)
                for lv = 1:len_inv
                    val = val+1;
                    for i = st:length(files)
                        D{i}.inv{val} = D{i}.inv{lv}; % duplicate
                        D{i}.inv{val}.inverse.trials = S.cond{cn};% specify
                    end
                end
            end

            % transform data
            if S.arcsinh_data
                for i = st:length(files)
                    files{i} = strrep(files{i},'.mat','_arcsinh.mat');
                    Dtrans{i} = clone(D{i}, files{i}, [D{i}.nchannels D{i}.nsamples D{i}.ntrials]);
                    Dtrans{i}(:) = log(D{i}(:)+sqrt(D{i}(:).^2+1));
                    D{i} = Dtrans{i};
                end
            end

            % Invert the forward model
            %==========================================================================
            if ~S.group_inversion
                for iv = 1:length(D{1}.inv)
                    for i = st:length(files)
                        D{i}.val=iv;
                        D{i} = spm_eeg_invert(D{i});
                    end
                end
            else
                for iv = 1:length(D{1}.inv)
                    for i = st:length(files)
                        D{i}.val=iv;
                    end
                    D = spm_eeg_invert(D);
                end
            end
            if ~iscell(D), D = {D}; end
        end



       % Save
        %==========================================================================
        for i = st:length(files)
            save(D{i});
        end
        

        % PRODUCE OUTPUTS
        if strcmp(S.sourceprior,'ECD')
            for i = st:length(files)
                % extract data
                xloc=D{i}.inv{end}.inverse.mniloc;
                xj=D{i}.inv{end}.inverse.jmni;
                for nd = 1:size(xloc{1},2) % n dipoles
                    temp=cat(3,xloc{:});
                    mniloc{g}(i,:,:,nd) = squeeze(temp(:,nd,:))';
                    temp=cat(3,xj{:});
                    jmni{g}(i,:,:,nd) = squeeze(temp(:,nd,:))';
                end
                lme{g}(i,:) = D{i}.inv{end}.inverse.F;
                pov{g}(i,:) = D{i}.inv{end}.inverse.pov;
            end
        else
            clear D
            % Compute conditional expectation of contrast and produce image
            %==========================================================================
            %mkdir('new_images')
            if ~isempty(images_out)
                for i = st:length(files)

                    D     = spm_eeg_load(fullfile(S.filepath,files{i}));

                    if isempty(S.imout)
                        imout = 1:size(images_out,1);
                    else
                        imout = S.imout;
                    end
                    %if max(imout) > 1%max(run_con)
                    %    for ii = imout(2:max(imout))
                    %        D.inv{ii} = D.inv{max(run_con)};
                    %    end
                    %elseif max(imout) > 1 && max(run_con)==1
                    %    for ii = imout
                    %        D.inv{ii} = D.inv{1};
                    %    end
                    %end

                    for con = imout

                        for sc = 1:length(S.cond)

                            woi=[];

                            % if it a baseline image?
                            if ischar(images_out{con,1})
                                if strcmp(images_out{con,1},'base')
                                    images_out{con,1} = basewin;
                                    basecon=con;
                                end
                            end

                            % specify D.val
                            if S.baseest 
                                if con==basecon
                                    D.val = (sc-1)*2+2;
                                else
                                    D.val = (sc-1)*2+1;
                                end
                            else
                                D.val = sc;
                            end
                            disp(['D_val = ' num2str(D.val)])

                            contrast = images_out{con,1};
                            if  size(contrast,2)==1   
                                if ~isempty(images_out{con,2})
                                    conwin = images_out{con,2}/2;
                                    if size(contrast,1)==1
                                        woi              = [contrast-conwin contrast+conwin];
                                    else 
                                        woi              = [contrast(grpind(i),1)-conwin contrast(grpind(i),1)+conwin];
                                    end
                                else
                                    woi              = contrast;
                                end
                            elseif size(contrast,2)==2 && size(contrast,1)>1
                                woi              = [contrast(grpind(i),1) contrast(grpind(i),2)];
                            else
                                woi              = contrast;
                            end

                            woi              = sort(woi)

                            %if max(woi)>max(contrasts{min(max(run_con),con),1}) || min(woi)<min(contrasts{min(max(run_con),con),1})
                            %    continue
                            %end

                            clear contrast
                            contrast.woi     = round([woi(1) woi(end)]);
                            %fboi             = 0;
                            %fboi             = sort(fboi);
                            %contrast.fboi    = round([fboi(1) fboi(end)]);
                            contrast.fboi    = S.freqwin;
                            contrast.display = 0;
                            contrast.smoothing  = S.smooth;
                            contrast.type = S.type;
                            D.inv{D.val}.contrast = contrast;

                            %cd(fullfile(swd,'SPM image files'))
                            D     = spm_eeg_inv_results(D);
                            D     = spm_eeg_inv_Mesh2Voxels_CAB(D); % CAB version saves by condition label (lines 204, 208 in function)
                        end
                    end
                end
            end

            % Cleanup
            %==========================================================================
            outfiles = dir(fullfile(S.filepath,'*spm*nii'));
            for f = 1:length(outfiles)
                movefile(fullfile(S.filepath,outfiles(f).name),S.outpath);

                if S.log_images
                    fname = fullfile(S.outpath,outfiles(f).name); 
                    spm_imcalc_ui(fname,fname,'log(i1+1)');
                end
            end
        end
    end % group
end % timewin
%end

if strcmp(S.sourceprior,'ECD')
    % save
    sname = ['ECD_out_' datestr(now,30)];
    save(fullfile(S.outpath,sname), 'mniloc','jmni','SubInd_ana','pov','lme','S');
end
