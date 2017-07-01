clear all
subjects = {'M1';'M2';'M4';'M6';'M7';'M10';'M12';'M14';'M15';'M16';'M17';'M29';'M30';'M32';'M33';'M9';'M19';'M20';'M21';'M22';'M24';'M25';'M26';'M28';'M35';'M36';'M37';'M38';'M40'};


fnames={'spm_base_';
  };


Ns = length(subjects);

run_con = [1];
contrasts = {[0 2000];
};
          
swd   = pwd;
    
for i = 1:Ns
    
subject = subjects(i);
subject = char(subject);   

fname=char(fnames);
fname= [fname subject '.mat'];

end

% Load data and set method
%==========================================================================
for i = 1:Ns
    
    fprintf('checking for previous inversions: subject %i\n',i);
    
    subject = subjects(i);
    subject = char(subject);   
    fname=char(fnames);
    fname= [fname subject '.mat'];
    
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
for i = 1:Ns
    fprintf('checking for foward models: subject %i\n',i);
    try
        [L, D{i}] = spm_eeg_lgainmat(D{i});
        Nd(i) = size(L,2);               % number of dipoles
    catch
        Nd(i) = 0;
    end
end

% use template head model where necessary
%==========================================================================
if max(Nd > 1024)
    NS = find(Nd ~= max(Nd));            % subjects requiring forward model
else
    NS = 1:Ns;
end
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

for i = 2:Ns
    [mod, list] = modality(D{i}, 1, 1);
    if ~ismember(inverse.modality, list)
        error([inverse.modality ' modality is missing from ' D{i}.fname]);
    end
end

% and save them (assume trials = types)
%--------------------------------------------------------------------------
for i = 1:Ns
    D{i}.inv{1}.inverse = inverse;
end

for i = NS

    fprintf('Registering and computing forward model (subject: %i)\n',i);

    % Forward model
    %----------------------------------------------------------------------
    D{i} = spm_eeg_inv_datareg_ui(D{i}, 1);
    D{i} = spm_eeg_inv_forward_ui(D{i}, 1);
    
    % save forward model
    %----------------------------------------------------------------------
    save(D{i});

end



for con = run_con
    
  for i = 1:Ns
      
    D{i}.val = con;
    D{i}.inv{con} = D{i}.inv{1}; 

%--------------------------------------------------------------------------
D{i}.inv{con}.inverse.trials = D{i}.condlist; % Trials
D{i}.inv{con}.inverse.type   = 'MSP';      % Priors on sources MSP, LOR or IID
D{i}.inv{con}.inverse.smooth = 1;        % Smoothness of source priors (mm)
D{i}.inv{con}.inverse.Np     = 256;         % Number of sparse priors (x 1/2)
D{i}.inv{con}.inverse.Han = 0;
D{i}.inv{con}.inverse.hpf = 250;
D{i}.inv{con}.inverse.lpf = 0;

contrast = contrasts(con);
contrast = [contrast{:}];
    if length(contrast) > 2
    woi              = [mean(contrast(i,1:2),2)-249 mean(contrast(i,1:2),2)+250];
    else
    woi              = contrast;
    end
    woi              = sort(woi);
    contrast.woi     = round([woi(1) woi(end)]);

D{i}.inv{con}.inverse.woi =woi; %1000*[min(D.time) max(D.time)];





    % get frequency window
    %----------------------------------------------------------------------
    fboi             = [90 120];
    fboi             = sort(fboi);
    contrast.fboi    = round([fboi(1) fboi(end)]);
    contrast.display = 0;
    contrast.smooth  = 10;
    
    str  = 'induced';

  end

  % Invert the forward model
%==========================================================================
D     = spm_eeg_invert(D);
if ~iscell(D), D = {D}; end

  
end



% Save
%==========================================================================
for i = 1:Ns
    save(D{i});
end
clear D

run_con = [1];
contrasts = {[0 2000];


% Compute conditional expectation of contrast and produce image
%==========================================================================
if ~isempty(contrast)

    % evaluate contrast and write image
    %----------------------------------------------------------------------
    for i = 1:Ns
        subject = subjects(i);
        subject = char(subject);   
        fname=char(fnames);
        fname= [fname subject '.mat'];
        D     = spm_eeg_load(fullfile(pwd,fname));
        
        for con = [1]
        D.val = con;
        contrast = contrasts(con);
        contrast = [contrast{:}];
        D.inv{con}.contrast = contrast;
        if length(contrast) > 2
            woi              = [contrast(i,1)-40 contrast(i,1)+40];
        else
            woi              = contrast;
        end
            woi              = sort(woi);
            contrast.woi     = round([woi(1) woi(end)]);
            D.inv{con}.contrast.woi = contrast.woi;
            D.inv{con}.inverse.woi = contrast.woi;
        D     = spm_eeg_inv_results(D);
        D     = spm_eeg_inv_Mesh2Voxels(D);
        end
    end
end

% Cleanup
%==========================================================================
cd(swd);
