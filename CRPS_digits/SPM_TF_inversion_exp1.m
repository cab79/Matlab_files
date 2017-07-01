clear all;

filepath = 'C:\Data\CRPS-DP\CRPS_Digit_Perception_exp1\correcttrials\';
run('M:\Matlab\Matlab_files\CRPS_digits\loadsubj.m');

cd(filepath);
%grplists = {33; 34; 31; 32} %unflipped,
%grplists = {[1];[29];[2];[30]}; %unflipped, left and right stimulation in different groups
grplists = {35; 36; 37; 38}; % Exp1 unflipped, left and right stimulation in different groups
%grplists = {39; 40; 41; 42}; %sublist_side = {'L','R','L','R'}; %Affected vs unaffected exp1
%lat = load('latencies_amp.dat','ACSII');
%ran = load('ranges_gfp.dat','ACSII');
add_baseline = 200;

ngrps = 4;
subspergrp = 13;



run_con = [1];
contrasts = {
    [-200 800],[];
            %[-200 0],[];
};

load freqs_limits
freq_peak_select = {
                     1, [1 2 3];
                     2, [1 2];   
                     };
%imout = [1 6 8];
%images_out = {
%       [-200 0],[];
%       [-200 0],[];
%    40,0;
%    88,10;
%    88,0;
%    [124 132],[];
%    160,0;
%    268,0;
%    144,0;
%    [172 324],[];
%    [544 600],[];
%};




templist = 1:subspergrp;

for g = 1:ngrps
    clear functions D;
    
    grpind = (g-1)*subspergrp+templist

    grplist = [grplists{g,:}]

    subjects_lists = subjlists(grplist);
    
    subjects = cell(1,1);
    
    for gl = 1:length(subjects_lists)
        subjects{1,1} = [subjects{1,1}; subjects_lists{gl,1}]
    end



    cd(filepath)

    Ns=0;
    for s = 1:length(subjects)
        for s2 = 1:length(subjects{s,1}) 
            Ns=Ns+1;
            tmp_nme = subjects{s,1}{s2,1};
            tmp_nme = strrep(tmp_nme, '.left', '_left');
            tmp_nme = strrep(tmp_nme, '.Left', '_left');
            tmp_nme = strrep(tmp_nme, '.right', '_right');
            tmp_nme = strrep(tmp_nme, '.Right', '_right');
            tmp_nme = strrep(tmp_nme, '.flip', '_flip');
            tmp_nme = strrep(tmp_nme, '.Flip', '_flip');
            tmp_nme = strrep(tmp_nme, '.aff', '_aff');
            tmp_nme = strrep(tmp_nme, '.Aff', '_aff');
            tmp_nme = strrep(tmp_nme, '.Unaff', '_unaff');
            tmp_nme = strrep(tmp_nme, '.unaff', '_unaff');
            tmp_nme = strrep(tmp_nme, '_Left', '_left');
            tmp_nme = strrep(tmp_nme, '_Right', '_right');
            tmp_nme = strrep(tmp_nme, '_Flip', '_flip');
            tmp_nme = strrep(tmp_nme, '_Aff', '_aff');
            tmp_nme = strrep(tmp_nme, '_Unaff', '_unaff');
            tmp_nme = strrep(tmp_nme, '.Exp1', '_Exp1');
            if strfind(tmp_nme,'left')
                fnames{Ns} = ['spm8_' tmp_nme];
            else
                
                fnames{Ns} = ['spm8_' tmp_nme];
            end
        end
    end


    swd   = pwd;

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
    if max(Nd > 1024)
        NS = find(Nd ~= max(Nd));            % subjects requiring forward model
    else
        NS = 1:Ns;
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
            %load meegfid;
            %load newmrifid;
            useheadshape = 1;
            D{i} = spm_eeg_inv_datareg_ui(D{i},1, [],[], useheadshape);
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

        D{i}.val = con;
        D{i}.inv{con} = D{i}.inv{1}; 
        
        Ds=D{i};
        Ds=struct(Ds);
        trials = cellfun(@str2num, unique({Ds.trials(:).label}));
        if ~isempty(intersect(trials,[1:5]))
            trials = cellfun(@num2str, num2cell(sort(trials)), 'UniformOutput', false);
        elseif ~isempty(intersect(trials,[6:10]))
            trials = cellfun(@num2str, num2cell(sort(trials,'descend')), 'UniformOutput', false);
        end
                
        D{i}.inv{con}.inverse.trials = trials; % Trials
        D{i}.inv{con}.inverse.type   = 'MSP';      % Priors on sources MSP, LOR or IID
        %D{i}.inv{con}.inverse.smooth = 0;        % Smoothness of source priors (mm)
        D{i}.inv{con}.inverse.Np     = 256;         % Number of sparse priors (x 1/2)
        %D{i}.inv{con}.inverse.Han = 1; 
        %D{i}.inv{con}.inverse.hpf = 48; 
        %D{i}.inv{con}.inverse.lpf = 0; % 0 for anticipation, 1 for LEPs (unless already filtered).

        invwin = contrasts{con,1};
        if  size(invwin,2)==1   
            if ~isempty(contrasts{con,2})
                conwin = contrasts{con,2}/2;
                if size(invwin,1)==1
                    woi              = [invwin-conwin invwin+conwin];
                else 
                    woi              = [invwin(grpind(i),1)-conwin invwin(grpind(i),1)+conwin];
                end
            else 
                woi              = invwin;
            end
        elseif size(invwin,2)==2 && size(invwin,1)>1
            woi              = [invwin(grpind(i),1) invwin(grpind(i),2)];
        else
            woi              = invwin;
        end
        woi              = sort(round(woi)) + add_baseline;
      
        D{i}.inv{con}.inverse.woi =woi; %1000*[min(D.time) max(D.time)];

      end

      % Invert the forward model for single-subject analysis
    %==========================================================================
  %  for i = st:Ns
  %    D{i}     = spm_eeg_invert(D{i});
  %      if ~iscell(D), D = {D}; end
  %  end
    
    % Invert the forward model for group analysis
    %==========================================================================
      D     = spm_eeg_invert(D);
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
        %mkdir('new_images')
        
            for i = st:Ns

                fname = fnames{i};
                D     = spm_eeg_load(fullfile(swd,fname));

                %if max(imout) > 2 && max(run_con)==2
                %    D.inv{max(imout)} = D.inv{2};
                %    for ii = imout(1:end-1)
                %        D.inv{ii} = D.inv{1};
                %    end
                %elseif max(imout) > 1 && max(run_con)==1
                %    for ii = imout
                %        D.inv{ii} = D.inv{1};
                %    end
                %end

                for f = 1:size(freq_peak_select,1)
                    %D.val = con
                    
                    %contrast = images_out{con,1};
                    %if  size(contrast,2)==1   
                    %    if ~isempty(images_out{con,2})
                    %        conwin = images_out{con,2}/2;
                    %        if size(contrast,1)==1
                    %            woi              = [contrast-conwin contrast+conwin];
                    %        else 
                    %            woi              = [contrast(grpind(i),1)-conwin contrast(grpind(i),1)+conwin];
                    %        end
                    %    else
                    %        woi              = contrast;
                    %    end
                    %elseif size(contrast,2)==2 && size(contrast,1)>1
                    %    woi              = [contrast(grpind(i),1) contrast(grpind(i),2)];
                    %else
                    %    woi              = contrast;
                    %end
                    freq = freqs_limits{freq_peak_select{f,1},1};
                    allwoi = freqs_limits{freq_peak_select{f,1},2};
                    selectwoi = allwoi(:,freq_peak_select{f,2});
                    
                    for w = 1:size(selectwoi,2)
                        woi              = sort(selectwoi(:,w)') + add_baseline
                        contrast.woi     = round([woi(1) woi(end)]);
                        fboi             = freq;
                        fboi             = sort(fboi);
                        contrast.fboi    = [fboi(1) fboi(end)];
                        %contrast.fboi    = [];
                        contrast.display = 0;
                        contrast.smoothing  = 12;
                        contrast.type = 'induced';
                        D.inv{con}.contrast = contrast;

                        %cd(fullfile(swd,'SPM image files'))
                        D     = spm_eeg_inv_results(D);
                        D     = spm_eeg_inv_Mesh2Voxels(D);
                    end
                end
            end

        % Cleanup
        %==========================================================================
        cd(swd);
end
%Subtract_Mean_Imcalc
matlabmail
