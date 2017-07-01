clear all;
%grplists = {39; 40; 41; 42}; %sublist_side = {'L','R','L','R'}; %Affected vs unaffected exp1
grplists = {1;2;29;30}; 
ngrps = 4;
subspergrp = 13;

if isunix
    filepath = '/scratch/cb802/Data/CRPS_Digit_Perception';
    run('/scratch/cb802/Matlab_files/CRPS_digits/loadsubj.m');
else
    filepath = 'W:\Data\CRPS_Digit_Perception';
    run('W:\Matlab_files\CRPS_digits\loadsubj.m');
end

templist = 1:subspergrp;

for g = 1:ngrps
    clear functions D;
    
    grpind = (g-1)*subspergrp+templist

    grplist = grplists{g,:}

    subjects = subjlists(grplist);

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
            fnames{Ns} = ['spm8_' tmp_nme];
        end
    end


    swd   = pwd;

    % Artefact rejection
    %==========================================================================
    %load batch_artefact;
    %S = matlabbatch{1,1}.spm.meeg.preproc;
%    clear S;
%    S.badchanthresh = 0.2;
%    S.methods(1).channels = 'all';
%    S.methods(1).fun = 'jump';
%    S.methods(1).settings.threshold = 75;
%    S.methods(2).channels = 'all';
%    S.methods(2).fun = 'peak2peak';
%    S.methods(2).settings.threshold = 200;
    
    %S.artefact.External_list = 0;
    %S.artefact.out_list = [];
    %S.artefact.in_list = [];
    %S.artefact.Weighted = 0;
    %S.artefact.Check_Threshold = 0;
    
 %   for i = 1:Ns
 %       fprintf('Artefact rejection: subject %i\n',i);
 %       fname = fnames{i};
 %       S.D = fullfile(pwd,fname);
 %       D = spm_eeg_artefact(S);
 %       %D     = spm_eeg_load(fullfile(pwd,fname));
 %       %save(D);
 %   end
    
     % Convert to 3D images (2 space and 1 time dimension) 
    %==========================================================================
    clear S D;
    %S.images.fmt = 'channels';
    %S.images.elecs = 
    %S.images.region_no = 
    %S.images.freqs =  
    %S.images.t_win = 
    S.n = 32;
    S.interpolate_bad = 1;
    
    for i = 1:Ns
        fprintf('Converting to images: subject %i\n',i);
        fname = ['a' fnames{i}];
        S.D = fullfile(pwd,fname);
        S.mode = 'scalp x time';
        spm_eeg_convert2images(S)
    end
   
    
     % Calculate variance images and smooth
    %==========================================================================    
    %load batch_smooth
    %for i = 1:Ns
    %    fname = ['a' fnames{i}];
    %    
    %    if strfind(fname, 'left')
    %        trials = [1:5];
    %    elseif strfind(fname, 'right')
    %        trials = [6:10];
    %    end
    %        
    %    for j = 1:length(trials)
    %        clear fnmes
    %        nme = dir(fullfile(filepath,fname,['type_' num2str(trials(j))],'trial*.img'));
    %        for k = 1:length(nme)
    %            fnmes{k,1} = fullfile(filepath,fname,['type_' num2str(trials(j))],nme(k).name);
    %        end
    %        outname = fullfile(filepath,fname,['type_' num2str(trials(j))],'var.img');
    %        if ~exist(outname, 'file')
    %            Output = spm_imcalc_ui(fnmes,outname,'std(X)',{1,[],[],[]})
    %            matlabbatch{1,1}.spm.spatial.smooth.data = {fullfile(filepath,fname,['type_' num2str(trials(j))],'var.img')};
    %            spm_jobman('run_nogui',matlabbatch);
    %        end
    %    end
    % end
    
    
end
