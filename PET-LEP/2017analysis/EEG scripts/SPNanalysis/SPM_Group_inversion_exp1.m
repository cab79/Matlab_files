clear all;
filepath = 'C:\Data\PET-LEP\Preprocessed';
anapath = 'C:\Data\PET-LEP\SPM_source_images';
run('C:\Data\Matlab\Matlab_files\PET-LEP\EEG scripts\loadsubj.m');
cd(filepath);
grplists = {1,2,3,... % healthy
            4,5,6};   % patient
add_baseline = {[-3500 -3000], [-500 0]};
suffix='.mat';
ngrps = size(grplists,1);
%subspergrp = 13;

%define latencies
load('timefreq_limits_ERP_evoked_1_2');
latency =  timefreq_limits.limits_all;%
peakdef =  timefreq_limits.bins;%            % defines which peak the latencies refer to. MUST BE CELL ARRAY
frequency = 0;%4:2:20;%[10]; % empty if ERP, performs freq analysis if specified

run_con = 1:(length(peakdef{:})+length(add_baseline));
%contrasts = {
%    [100 400],[];
%};

contrasts = cell(1,1);
for b = 1:length(add_baseline)
    contrasts{b,1} = add_baseline{b};
    contrasts{b,2} = []; 
end
for c = 1:length(peakdef{:})
    contrasts{c+length(add_baseline),1} = [latency{1,1}{1,c}(1), latency{1,1}{1,c}(end)]*2-3502;
    contrasts{c+length(add_baseline),2} = []; 
end

baselines = cell(1,1);
for b = 1:length(add_baseline)
    baselines{b,1} = contrasts{b,1};
end
for b = 1:length(peakdef{:})
    if contrasts{b+length(add_baseline),1}(1)>0
        baselines{b+length(add_baseline),1} = [-500 0];
    else
        baselines{b+length(add_baseline),1} = [-3500 -3000];
    end
end

%imout = [1:2];
%images_out = {
%       [-200 0],[];
%       [-200 0],[];
%};

imout = run_con;
images_out = contrasts;

%templist = 1:subspergrp;
sourcepriors = 'LOR';
grpind = 'ind';

newdir = [sourcepriors '_' grpind '_' datestr(datetime,30)];

save(fullfile(anapath,[newdir 'imageinfo.mat']));

for g = 1:ngrps
    clear functions D;
    %grpind = (g-1)*subspergrp+templist
    grplist = [grplists{g,:}]
    subjects_lists = subjlists(grplist);
    subjects = cell(1,1);
    for gl = 1:size(subjects_lists,1)
        subjects{1,1} = [subjects{1,1}; subjects_lists{gl,1}]
    end

    cd(filepath)
    Ns=0;
    for s = 1:length(subjects)
        for s2 = 1:length(subjects{s,1}) 
            Ns=Ns+1;
            tmp_nme = subjects{s,1}{s2,1};
            fnames{Ns} = ['spm12_' tmp_nme suffix];
        end
    end

    % Load data and set method
    %==========================================================================
    %if g==2; st = 11; else st =1;end
    st=1;
    for i = st:Ns

        fprintf('checking for previous inversions: subject %i\n',i);

        fname = fnames{i};

        D{i}     = spm_eeg_load(fullfile(pwd,fname));
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
    Nd = zeros(1,Ns);
    for i = st:Ns
        fprintf('checking for foward models: subject %i\n',i);
        try
            [L, D{i}] = spm_eeg_lgainmat(D{i});
            Nd(i) = size(L,2);               % number of dipoles
        catch
            Nd(i) = 0;
        end
        clear L;
    end

    % use template head model where necessary
    %==========================================================================
    %if max(Nd > 1024)
    %    NS = find(Nd ~= max(Nd));            % subjects requiring forward model
    %else
        NS = 1:Ns;
    %end
    
    
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
            load meegfid_new2;
            load newmrifid;
            %D{i}.fiducials.pnt(2:3,1) = -D{i}.fiducials.pnt(2:3,1);
            useheadshape = 1;
            D{i} = spm_eeg_inv_datareg_ui(D{i},1, meegfid,newmrifid, useheadshape);
            %D{i} = spm_eeg_inv_datareg_ui(D{i},1);
            D{i} = spm_eeg_inv_forward_ui(D{i}, 1);
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

    for con = run_con

      for i = st:Ns
          
          % baseline correct
          if ~isempty(baselines{con,1})
              %S.D = D{i};
              %S.timewin = baselines{con,1};
              %S.save=0;
              %D{i} = spm_eeg_bc(S);
              t=baselines{con,1};
              t(1) = D{i}.indsample(baselines{con,1}(1)/1000);
              t(2) = D{i}.indsample(baselines{con,1}(2)/1000);
                indchannels = [D{i}.indchantype('Filtered') D{i}.indchantype('MEGCOMB')];
                for k = 1: D{i}.ntrials
                    tmp = mean(D{i}(indchannels, t(1):t(2), k), 2);
                    D{i}(indchannels, :, k) = D{i}(indchannels, :, k) - repmat(tmp, 1, D{i}.nsamples);
                end
          end
          

        D{i}.val = con;
        D{i}.inv{con} = D{i}.inv{1}; 
        D{i}.inv{con}.inverse.trials = D{i}.condlist; % Trials
        D{i}.inv{con}.inverse.type   = sourcepriors;      % Priors on sources MSP, LOR or IID
        %D{i}.inv{con}.inverse.smooth = 0;        % Smoothness of source priors (mm)
        D{i}.inv{con}.inverse.Np     = 256;         % Number of sparse priors (x 1/2)
        D{i}.inv{con}.inverse.Han = 1; 
        %D{i}.inv{con}.inverse.hpf = 48; 
        %D{i}.inv{con}.inverse.lpf = 0; % 0 for anticipation, 1 for LEPs (unless already filtered).

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
        if strcmp(grpind,'grp')
            D   = spm_eeg_invert(D);
        elseif strcmp(grpind,'ind')
            for i = st:Ns
                D{i}   = spm_eeg_invert(D{i});
            end
        end
                
        if ~iscell(D), D = {D}; end

    end



   % Save
%==========================================================================
    for i = st:Ns
        save(D{i});
    end
    clear D


    % Compute conditional expectation of contrast and produce image
    %==========================================================================

    if ~isempty(images_out)
        for i = st:Ns

            fname = fnames{i};
            D     = spm_eeg_load(fullfile(filepath,fname));

            if max(imout) > 2 && max(run_con)==2
                D.inv{max(imout)} = D.inv{2};
                for ii = imout(1:end-1)
                    D.inv{ii} = D.inv{1};
                end
            elseif max(imout) > 1 && max(run_con)==1
                for ii = imout
                    D.inv{ii} = D.inv{1};
                end
            end

            for con = imout
                D.val = con

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
                clear contrast
                contrast.woi     = round([woi(1) woi(end)]);
                %fboi             = 0;
                %fboi             = sort(fboi);
                %contrast.fboi    = round([fboi(1) fboi(end)]);
                contrast.fboi    = [];
                contrast.display = 0;
                contrast.smoothing  = 8;
                contrast.type = 'evoked';
                D.inv{con}.contrast = contrast;

                D     = spm_eeg_inv_results(D);
                D     = spm_eeg_inv_Mesh2Voxels(D);
            end
        end
    end
    
    if ~exist(fullfile(anapath,newdir),'dir')
        mkdir(fullfile(anapath,newdir));
    end
    movefile(fullfile(filepath,'*.nii'),fullfile(anapath,newdir));

    % Cleanup
    %==========================================================================
    cd(filepath);
end
%Subtract_Mean_Imcalc
matlabmail
