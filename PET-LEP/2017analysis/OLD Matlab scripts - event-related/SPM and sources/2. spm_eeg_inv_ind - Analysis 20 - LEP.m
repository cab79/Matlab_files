% CURRENT PERIOD COVERS n2 AND p2

clear all
subjects = {'P1_';'P2_';'P5_';'P6_';'P7_';'P8_';'P14_';'P16_';'P20_';'P22_';'P23_';'P24_';'P25_';'P27_';'P30_';'P31_';'P32_';'P33_';'P35_';'S1_';'S2_';'S3_';'S5_';'S6_';'S8_';'S9_';'S10_';'S11_';'S20_'};



sub_ana = 1:length(subjects);
%sub_ana = 1:15;

fnames={'mspm_';
  };


load N1latencies
load N2latencies
load P2latencies

Ns = length(subjects);


          
swd   = pwd;
    
%for i = sub_ana
%subject = subjects(i);
%subject = char(subject);   
%fname=char(fnames);
%fname= [fname subject '.mat'];
%end

% Load data and set method
%==========================================================================
for i = sub_ana
    run_con = [6];
    contrasts = {[-3500 -3000];
        [-3500 0];
        [-3500 0];
        [-500 0];
        [-200 800];
        [-200 800];
    };
    
    fprintf('checking for previous inversions: subject %i\n',i);
    
    subject = subjects(i);
    subject = char(subject);   
    fname=char(fnames);
    fname= [fname subject '_LEP.mat'];
    
    D     = spm_eeg_load(fullfile(pwd,fname));
    S.D = D;
    S.time = [-200 0];
    S.save = 0;
    D = spm_eeg_bc(S);
    
    D.val             = 1;
    D.inv{1}.method   = 'Imaging';
    
    % clear redundant models
    %----------------------------------------------------------------------
    D.inv = D.inv(1);
    
    
    % clear previous inversions
    %----------------------------------------------------------------------
    try, D.inv{1} = rmfield(D.inv{1},'inverse' ); end
    try, D.inv{1} = rmfield(D.inv{1},'contrast'); end
    
    % save forward model parameters
    %----------------------------------------------------------------------
    save(D);
    
%end

% Check for existing forward models and consistent Gain matrices
%--------------------------------------------------------------------------
%Nd = zeros(1,length(sub_ana));
%for i = sub_ana
    fprintf('checking for foward models: subject %i\n',i);
    try
        [L, D] = spm_eeg_lgainmat(D);
        Nd = size(L,2);               % number of dipoles
    catch
        Nd = 0;
    end
%end

% use template head model where necessary
%==========================================================================
if Nd==0
%    NS = find(Nd ~= max(Nd));            % subjects requiring forward model
%else
%    NS = sub_ana;
%end

%for i = NS

    cd(D.path);

    % specify cortical mesh size (1 tp 4; 1 = 5125, 2 = 8196 dipoles)
    %----------------------------------------------------------------------
    Msize  = 1;

    % use a template head model and associated meshes
    %======================================================================
    D = spm_eeg_inv_mesh_ui(D, 1, 1, Msize);

    % save forward model parameters
    %----------------------------------------------------------------------
    save(D);

end

% Select modality
%==========================================================================
inverse.modality  = spm_eeg_modality_ui(D, 1, 1);

%for i = sub_ana(2:end)
    [mod, list] = modality(D, 1, 1);
    if ~ismember(inverse.modality, list)
        error([inverse.modality ' modality is missing from ' D.fname]);
    end
%end

% and save them (assume trials = types)
%--------------------------------------------------------------------------
%for i = sub_ana
    D.inv{1}.inverse = inverse;
%end

if Nd==0

    fprintf('Registering and computing forward model (subject: %i)\n',i);

    % Forward model
    %----------------------------------------------------------------------
    %D = spm_eeg_inv_datareg_ui(D, 1);
    load meegfid;
    load newmrifid;
    useheadshape = 0;
    D = spm_eeg_inv_datareg_ui(D,1, meegfid, newmrifid, useheadshape);
    %D = spm_eeg_inv_datareg_ui(D,1);
    D = spm_eeg_inv_forward_ui(D, 1);
    
    
    % save forward model
    %----------------------------------------------------------------------
    save(D);

end



for con = run_con
    
  %for i = sub_ana
      
    D.val = con;
    D.inv{con} = D.inv{1}; 
    D.inv{con}.inverse.trials = D.condlist; % Trials
    D.inv{con}.inverse.type   = 'LOR';      % Priors on sources MSP, LOR or IID
    D.inv{con}.inverse.smooth = 1.0;        % Smoothness of source priors (mm) % Use value of 1.
    D.inv{con}.inverse.Np     = 128;         % Number of sparse priors (x 1/2)
    D.inv{con}.inverse.Han = 0; % Leave at 0.
    D.inv{con}.inverse.hpf = 20; 
    D.inv{con}.inverse.lpf = 1; % 0 for anticipation, 1 for LEPs.

    contrast = contrasts(con);
    contrast = [contrast{:}];
    if length(contrast) > 2
    woi              = [mean(contrast(i,1:2),2)-249 mean(contrast(i,1:2),2)+250];
    else
    woi              = contrast;
    end
    woi              = sort(woi);
    %contrast.woi     = round([woi(1) woi(end)]);

    D.inv{con}.inverse.woi =woi; %1000*[min(D.time) max(D.time)];

  %end

  % Invert the forward model
%==========================================================================
  D     = spm_eeg_invert(D);
  %if ~iscell(D), D = {D}; end

  
end



% Save
%==========================================================================
%for i = sub_ana
    save(D);
%end
%clear D

run_con = [6];
contrasts = {[-3500 -3000];
    [-3500 -3000];
    [-500 0];
    [-500 0];
    [-200 0];
    [380 420];
    %[min([N2LAT(:,1) N2LAT(:,2)])-50 max([N2LAT(:,1) N2LAT(:,2)])+50];
    N1LAT(:,1);
    N1LAT(:,2);
    N2LAT(:,1);
    N2LAT(:,2);
    P2LAT(:,1);
    P2LAT(:,2);
    };

% Compute conditional expectation of contrast and produce image
%==========================================================================
if ~isempty(contrast)
    %for i = sub_ana
        
        
        
        subject = subjects(i);
        subject = char(subject);   
        fname=char(fnames);
        fname= [fname subject '_LEP.mat'];
        D     = spm_eeg_load(fullfile(pwd,fname));
        
        %for con = [7:11]
        %    D.inv{con} = D.inv{6};
        %end
        
        for con = run_con
            D.val = con;
            contrast = contrasts(con);
            contrast = [contrast{:}];
            %D.inv{con}.contrast = contrast;
            if length(contrast) > 2
                woi              = [contrast(i,1)-20 contrast(i,1)+20];
            else
                woi              = contrast;
            end
            woi              = sort(woi);
            contrast.woi     = round([woi(1) woi(end)]);
            fboi             = 0;
            fboi             = sort(fboi);
            contrast.fboi    = round([fboi(1) fboi(end)]);
            contrast.display = 0;
            contrast.smoothing  = 8;
            contrast.type = 'evoked';
            D.inv{con}.contrast = contrast;
            D     = spm_eeg_inv_results(D);
            D     = spm_eeg_inv_Mesh2Voxels(D);
        end
    %end
end
end
% Cleanup
%==========================================================================
cd(swd);
