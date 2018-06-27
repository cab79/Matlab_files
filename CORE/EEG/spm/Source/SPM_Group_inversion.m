%% SPM source analysis: setup

clear all;

%load('C:\Data\Catastrophising study\meegfid\fidmat2.mat');
%for fm = 1:length(fidmat)

%% generic directories for all analyses for this study
%-------------------------------------------------------------
% root directory in which SPM data files are located
S.filepath = 'C:\Data\Catastrophising study\SPMdata'; 
% place to save source images
S.outpath = 'C:\Data\Catastrophising study\SPMdata\sourceimages_GS'; 
% load .xlsx file containing 'Participant_ID', 'Group', and covariates
S.pdatfile = 'C:\Data\Catastrophising study\Behavioural\Participant_data_nocodes.xlsx';
%fiducials directory
S.fid_dir='C:\Data\Catastrophising study\meegfid';

%% specific directory and file information for this analysis
%-------------------------------------------------------------
% prefix, middle part, or suffix of files to load (or leave empty) to select a subset of files in
% the folder
S.fpref = 'mspm12';
S.fmid = '';
%S.fsuff = '_orig_cleaned.mat';
S.fsuff = '_orig_cleaned_trialNmatch.mat';

%% specific settings for this analysis
% which codes to analyse in 'Include' columns in participant data file?
S.include_codes = [1];
S.grps = {1;2}; %inversion on each group separately: separate with colon. Otherwise separate with comma
% time and frequecy windows
S.freqwin = []; % empty if not requiring freq analysis
%S.timewin = [-5500 -2500]; % Best to narrow this as much as possible to the range of interest. Empty will include whole epoch. 
%S.basewin = [-5500 -5000]; % empty will not baseline correct and will not produce a baseline image.
%S.timewin = [-3000 0]; % empty will include whole epoch. Best to narrow this as much as possible to the range of interest.
%S.basewin = [-3000 -2500]; % empty will not baseline correct and will not produce a baseline image.
S.timewin = [100 800]; % empty will include whole epoch. Best to narrow this as much as possible to the range of interest.%
S.basewin = [-500 0]; % empty will not baseline correct and will not produce a baseline image.
S.baseest = 1; % 1 = baseline estimated in separate model
% latency ranges of images to procedure. Will only produce images within
% the S.timewin range.
imout = []; % select which to produce on this run, or leave empty to run all
images_out = {
        'base',[]; % outputs a baseline image at the range set by S.basewin
   %     [-4660 -4278],[]; %Exp
   %     [-4564 -4506],[]; %Grp
   %     [-4552 -4468],[]; %Grp*Exp
   %     [-4454 -4376],[]; %Grp
        %[-4308 -4246],[]; %Grp*Exp
        %[-3364 -3292],[]; %Grp*Exp
        %[-3264 -3210],[]; %Grp*Exp
        %[-3050 -2946],[]; %Exp
   %     [-2698 -2500],[]; %Exp, Grp*Exp
        %[-2514 -2448],[]; %Exp
   %     [-2264 -2226],[]; %Exp
   %     [-2252 -2202],[]; %Grp
        %[-2162 -2002],[]; %Grp
        %[-1288 -1076],[]; %Exp
        [416 478],[]; %Exp
};

%smooth output images (specify FWHM or 0 for no smoothing)
S.smooth = 12;
S.type = 'evoked';
S.sourceprior = 'GS'; % Priors on sources, e.g. MSP, GS, LOR or IID
S.Npriors = 256; % Number of sparse priors (x 1/2 brain)
S.Han = 1; % apply Hanning window

%% RUN

%S.outpath = [S.outpath num2str(fm)];
if ~exist(S.outpath,'dir')
    mkdir(S.outpath);
end

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

for g = 1:Ngrp_ana
    subID={};
    files={};
    Nsub_ana(g) = length(SubInd_ana{g});
    for s = 1:Nsub_ana(g)
        subID{s} = pdata{SubInd_ana{g}(s)+1,sub_col};
        if isnumeric(subID{s}); subID{s} = num2str(subID{s}); end;
        fname = dir(fullfile(S.filepath,[S.fpref '*' subID{s} '*' S.fmid  '*' S.fsuff]));
        files{s} = fname.name;
        % Baseline Correction
        clear B D
        B.D = fullfile(S.filepath,files{s});
        B.timewin = S.basewin;
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
        
        D{i}.val             = 1;
        D{i}.inv{1}.method   = 'Imaging';

        % clear redundant models
        %----------------------------------------------------------------------
        D{i}.inv = D{i}.inv(1);


        % clear previous inversions
        %----------------------------------------------------------------------
        try, D{i}.inv{1} = rmfield(D{i}.inv{1},'inverse' ); end
        try, D{i}.inv{1} = rmfield(D{i}.inv{1},'contrast'); end

        % save forward model parameters
        %----------------------------------------------------------------------
        save(D{i});

    end
    
    % Check for existing forward models and consistent Gain matrices
    %--------------------------------------------------------------------------
    Nd = zeros(1,length(st:length(files)));
    for i = st:length(files)
        fprintf('checking for foward models: subject %i\n',i);
        %try
        %    [L, D{i}] = spm_eeg_lgainmat(D{i});
        %    Nd(i) = size(L,2);               % number of dipoles
        %catch
            Nd(i) = 0;
        %end
        clear L;
    end

    % use template head model where necessary
    %==========================================================================
    if max(Nd > 1024)
        NS = find(Nd ~= max(Nd));            % subjects requiring forward model
    else
        NS = st:length(files);
    end
    
    
    if ~isempty(NS)
        for i = NS

            cd(D{i}.path);

            % specify cortical mesh size (1 tp 4; 1 = 5125, 2 = 8196 dipoles)
            %----------------------------------------------------------------------
            Msize  = 2;

            % use a template head model and associated meshes
            %======================================================================
            D{i} = spm_eeg_inv_mesh_ui(D{i}, 1, 1, Msize);

            % save forward model parameters
            %----------------------------------------------------------------------
            save(D{i});
        end

        % Select modality
        %==========================================================================
        inverse.modality  = spm_eeg_modality_ui(D{1}, 1, 1);

        for i = NS(2:end)
            [mod, list] = modality(D{i}, 1, 1);
            if ~ismember(inverse.modality, list)
                error([inverse.modality ' modality is missing from ' D{i}.fname]);
            end
        end

        % and save them (assume trials = types)
        %--------------------------------------------------------------------------
        for i = NS
            D{i}.inv{1}.inverse = inverse;
        end


        for i = NS(1)

            %if ~exist('D','var'); load tempD; end;

            fprintf('Registering and computing forward model (subject: %i)\n',i);

            % Forward model
            %----------------------------------------------------------------------
            %D{i} = spm_eeg_inv_datareg_ui(D{i}, 1);
            load(fullfile(S.fid_dir,'meegfid2.mat'));
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
            D{i} = spm_eeg_inv_datareg_ui(D{i},1, meegfid,newmrifid, useheadshape);
            %D{i} = spm_eeg_inv_datareg_ui(D{i},1);
            D{i} = spm_eeg_inv_forward_ui_CAB(D{i}, 1);
            clear functions
            [L, D{i}] = spm_eeg_lgainmat(D{i});
            clear functions
            % save forward model
            %----------------------------------------------------------------------
            save(D{i});
            clear L;

        end

        for i = NS(2:end)

            D{i}.inv{1} = D{NS(1)}.inv{1};
            %D{i}.inv{1} = rmfield(D{i}.inv{1},'gainmat');

            %if ~exist('D','var'); load tempD; end;

            fprintf('Registering and computing forward model (subject: %i)\n',i);

            %clear functions
            [L, D{i}] = spm_eeg_lgainmat(D{i});
            %clear functions
            % save forward model
            %----------------------------------------------------------------------
            save(D{i});
            clear L;

        end
    end

    if ~exist('contrasts','var')
        if S.baseest
            contrasts = {S.basewin,[]; S.timewin,[]};
        else
            contrasts = {S.timewin,[]};
        end
    end
    
    if ~exist('run_con','var')
        if S.baseest
            run_con = 1:2;
        else
            run_con = 1;
        end
    end
    for con = run_con

      for i = st:length(files)

        D{i}.val = con;
        D{i}.inv{con} = D{i}.inv{1}; 
        D{i}.inv{con}.inverse.trials = D{i}.condlist; % Trials
        D{i}.inv{con}.inverse.type   = S.sourceprior;
        %D{i}.inv{con}.inverse.smooth = 0;        % Smoothness of source priors (mm)
        D{i}.inv{con}.inverse.Np     = S.Npriors;
        D{i}.inv{con}.inverse.Han = S.Han; 
        if ~isempty(S.freqwin)
            D{i}.inv{con}.inverse.hpf = S.freqwin(2); 
            D{i}.inv{con}.inverse.lpf = S.freqwin(1); 
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

        D{i}.inv{con}.inverse.woi =woi; %1000*[min(D.time) max(D.time)];

      end

      % Invert the forward model
    %==========================================================================

      D     = spm_eeg_invert(D);
        if ~iscell(D), D = {D}; end

    end



   % Save
    %==========================================================================
    for i = st:length(files)
        save(D{i});
    end
    clear D


    % Compute conditional expectation of contrast and produce image
    %==========================================================================
    %mkdir('new_images')
    if ~isempty(images_out)
        for i = st:length(files)

            D     = spm_eeg_load(fullfile(S.filepath,files{i}));

            if isempty(imout)
                imout = 1:size(images_out,1);
            end
            if max(imout) > max(run_con)
                for ii = imout(max(run_con)+1:max(imout))
                    D.inv{ii} = D.inv{max(run_con)};
                end
            %elseif max(imout) > 1 && max(run_con)==1
            %    for ii = imout
            %        D.inv{ii} = D.inv{1};
            %    end
            end

            for con = imout
                D.val = con

                if ischar(images_out{con,1})
                    if strcmp(images_out{con,1},'base')
                        images_out{con,1} = S.basewin;
                    end
                end
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
                
                if max(woi)>max(contrasts{con,1}) || min(woi)<min(contrasts{con,1})
                    continue
                end
                
                clear contrast
                contrast.woi     = round([woi(1) woi(end)]);
                %fboi             = 0;
                %fboi             = sort(fboi);
                %contrast.fboi    = round([fboi(1) fboi(end)]);
                contrast.fboi    = S.freqwin;
                contrast.display = 0;
                contrast.smoothing  = S.smooth;
                contrast.type = S.type;
                D.inv{con}.contrast = contrast;

                %cd(fullfile(swd,'SPM image files'))
                D     = spm_eeg_inv_results(D);
                D     = spm_eeg_inv_Mesh2Voxels_CAB(D); % CAB version saves by condition label (lines 204, 208 in function)
            end
        end
    end

    % Cleanup
    %==========================================================================
    outfiles = dir(fullfile(S.filepath,'mspm*nii'));
    for f = 1:length(outfiles)
        movefile(fullfile(S.filepath,outfiles(f).name),S.outpath);
    end
end
%end