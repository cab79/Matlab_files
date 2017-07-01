clear all;
%filepath = 'C:\Data\CRPS-DP\CRPS_Digit_Perception\SPM image files\eleposnorm';
filepath = 'C:\Data\CRPS-DP\CRPS_Digit_Perception_exp1\alltrials\SPM image files\eleposnorm';
run('C:\Data\Matlab\Matlab_files\CRPS_digits\loadsubj.m');
cd(filepath);
%grplists = {[1],[29];[2],[30]}; %unflipped, left and right stimulation in different groups
grplists = {35, 37; 36, 38}; % Exp1 unflipped, left and right stimulation in different groups
%lat = load('latencies_amp.dat','ACSII');
%ran = load('ranges_gfp.dat','ACSII');
baseline = [0 200];
add_baseline = 200;
suffix='';
use_flipped = 0; %SPM SOURCE ANALYSIS CANNOT INTERPRET THE FLIPPED CHANNEL POSITIONS
ngrps = 4;
subspergrp = 13;

%% ind lat
%run_con = [1:8];
%contrasts = {
%       [-200 0],[];
%    lat(:,1),20;
%    lat(:,2),40;
%    lat(:,3),40;
%    lat(:,4),40;
%    lat(:,5),40;
%    lat(:,6),40;
%    lat(:,7),40;
%};

%% mean lat
%run_con = [1:8];
%lat = mean(lat,1);
%contrasts = {
%       [-200 0],[];
%    lat(1),20;
%    lat(2),28;
%    lat(3),40;
%    lat(4),32;
%    lat(5),48;
%    lat(6),152;
%    lat(7),196;
%};
%run_con = [5:6];
%contrasts = {
%       [-200 0],[];
%    [0 200],[];
%    [200 400],[];
%    180,20;
%    [-200 400],[];
%    [-200 200],[];
%    [-200 200],[];
%};

%run_con = [2 6];
%contrasts = {
%       [-200 400],[];
%       [-200 0],[];
%       [0 400],[];
%       [0 400],[];
%       88,40;
%       [0 400],[];
%       [0 200],[];
%};

%imout = [2 6];
%images_out = {
%       [-200 0],[];
%       [-200 0],[];
%       [0 400],[];
%       40,10;
%       88,10;
%       88,10;
%       88,10;
%};

run_con = [1:3];
contrasts = {
    [-200 0],[];
    [-200 400],[];
    [100 400],[];
       %[0 200],[];
       %[0 100],[];
       %[0 100],[];
       %[-200 400],[];
       %[0 200],[];
       %[0 200],[];
       %[0 400],[];
};

imout = [1:4];
images_out = {
       [-200 0],[];
       [-200 0],[];
    132,0;
    [124 140],[];
    [120 180],[];
    [124 132],[];
    [216 324],[];
};

%run_con = [1:2];
%contrasts = {
%    [-200 800],[];
%      [-200 0],[];
%};


%imout = [1:8];
%images_out = {
%    [-200 0]
%    [0 20];
%    [20 70];
%    [70 106];
%    [106 144];
%    [144 360];
%    [400 800];
%};

%run_con = [4:6];
%contrasts = {
%        lat(:,1),20;
%    lat(:,2),40;
%    lat(:,3),40;
%    lat(:,4),40;
%    lat(:,5),40;
%    lat(:,6),40;
%    lat(:,7),40;
%       [-200 0],[];
%};

%imout = [4 5 6 8];
%images_out = {
%        lat(:,1),20;
%    lat(:,2),40;
%    lat(:,3),40;
%    lat(:,4),40;
%    lat(:,5),40;
%    lat(:,6),40;
%    lat(:,7),40;
%       [-200 0],[];
%};

%imout = [4 5 6 8];
%images_out = {
%        ran(:,1:2);
%    ran(:,3:4);
%    ran(:,5:6);
%    ran(:,7:8);
%    ran(:,9:10);
%    ran(:,11:12);
%    ran(:,13:14);
%       [-200 0];
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
            if use_flipped==1
                if strfind(tmp_nme, 'left')
                    fnames{Ns} = ['maspm8_' tmp_nme suffix];
                elseif strfind(tmp_nme, 'right')
                    fnames{Ns} = ['maspm8_flip_' tmp_nme suffix];
                end
            else
                if strfind(tmp_nme, 'left')
                elseif strfind(tmp_nme, 'right')
                end
                fnames{Ns} = ['maspm8_' tmp_nme suffix];
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
        
        S=struct;
        S.D=D{i};
        S.timewin = baseline;
        S.save = 0;
        D{i} = spm_eeg_bc(S);
        
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
            %D{i}.fiducials.pnt(2:3,1) = -D{i}.fiducials.pnt(2:3,1);
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
        D{i}.inv{con}.inverse.trials = D{i}.condlist; % Trials
        D{i}.inv{con}.inverse.type   = 'GS';      % Priors on sources MSP, LOR or IID
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
        woi              = sort(woi) + add_baseline;
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
for i = st:Ns
    save(D{i});
end
clear D


        % Compute conditional expectation of contrast and produce image
        %==========================================================================
        %mkdir('new_images')
        if ~isempty(images_out)
            for i = st:Ns

                fname = fnames{i};
                D     = spm_eeg_load(fullfile(swd,fname));

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
                    
                    woi              = sort(woi) + add_baseline
                    clear contrast
                    contrast.woi     = round([woi(1) woi(end)]);
                    %fboi             = 0;
                    %fboi             = sort(fboi);
                    %contrast.fboi    = round([fboi(1) fboi(end)]);
                    contrast.fboi    = [];
                    contrast.display = 0;
                    contrast.smoothing  = 12;
                    contrast.type = 'evoked';
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
