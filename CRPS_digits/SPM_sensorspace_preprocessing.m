clear all;

   %filepath = 'C:\Data\CRPS-DP\CRPS_Digit_Perception\SPM image files\eleposnorm';
   filepath = 'C:\Data\CRPS-DP\CRPS_Digit_Perception_exp1\alltrials';
    run('C:\Data\Matlab\Matlab_files\CRPS_digits\loadsubj.m');
%end
cd(filepath);
%grplists = {36;38}; %sublist_side = {'L','R','L','R'}; %Affected vs unaffected exp1
grplists = {39; 40; 41; 42}; %sublist_side = {'L','R','L','R'}; %Affected vs unaffected exp1
%grplists = {2; 30}; %sublist_side = {'L','R','L','R'}; %Affected vs unaffected exp2
%grplists = {47;48;49;50}; %sublist_side = {'L','R','L','R'}; %Affected vs unaffected exp2
%grplists = {35;37};
%grplists = {37};
ngrps = length(grplists);
subspergrp = 13;
templist = 1:subspergrp;
start = 1;
use_flipped=1;
balance_trials = 1;
suffix = '';
%suffix = '_change';

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
            tmp_nme = strrep(tmp_nme, '.Exp2', '_Exp2');
            
            %if strfind(tmp_nme, 'left'); continue;end;
            
            if use_flipped==1
                if strfind(tmp_nme, 'left')
                    fnames{Ns} = ['spm8_' tmp_nme suffix];
                elseif strfind(tmp_nme, 'right')
                    fnames{Ns} = ['spm8_flip_' tmp_nme suffix];
                end
            else
                fnames{Ns} = ['spm8_' tmp_nme suffix];
            end
            
            
        end
    end
    swd   = pwd;
    for i = 1:subspergrp
        % Artefact rejection
        %==========================================================================
        %load batch_artefact;
        %S = matlabbatch{1,1}.spm.meeg.preproc;
        clear S;
        S.badchanthresh = 0.2;
        S.methods(1).channels = 'all';
        S.methods(1).fun = 'jump';
        S.methods(1).settings.threshold = 75;
        S.methods(2).channels = 'all';
        S.methods(2).fun = 'peak2peak';
        S.methods(2).settings.threshold = 200;
        %S.artefact.External_list = 0;
        %S.artefact.out_list = [];
        %S.artefact.in_list = [];
        %S.artefact.Weighted = 0;
        %S.artefact.Check_Threshold = 0;
       % for i = 1:subspergrp
            fprintf('Artefact rejection: subject %i\n',i);
            fname = fnames{i};
            %if (~isempty(strfind(fname, 'right')) || ~isempty(strfind(fname, 'Right'))); 
            %    continue
            %end
            S.D = fullfile(pwd,fname)
            D = spm_eeg_artefact(S);
            %D     = spm_eeg_load(fullfile(pwd,fname));
            %save(D);
       % end%%
       
       
        % Match imbalanced conditions
        %==========================================================================
        if balance_trials
            balance = {{'1','2','3','4','5'; % old conditions
                        1 2 2 2 1}; % new conditions to balance with respect to - number of trials of 1 will be the same as ntrials of 2
                        {'6','7','8','9','10';
                        1 2 2 2 1};
                        };
            E=struct(D);
            [A,B,C] = unique({E.trials.label});
            for u = 1:length(A)
                nA(u) = sum(C==u);
            end
            % get new conditions to balance with
            sumtypes = [];
            if any(ismember(A,balance{1}(1,:)))
                for ii = 1:length(A)
                    sumtypes(ii) = balance{1}{2,strcmp(A{ii},balance{1}(1,:))};
                end
            elseif any(ismember(A,balance{2}(1,:)))
                for ii = 1:length(A)
                    sumtypes(ii) = balance{2}{2,strcmp(A{ii},balance{2}(1,:))};
                end
            end
            nC = unique(sumtypes);
            for ii = nC
                nT(ii) = sum(nA(sumtypes==ii));
            end
            maxi = find(ismember(C,find(sumtypes==nC(find(nT==max(nT)))))); % indices of C that are from the conditions with max number of trials
            mini = find(ismember(C,find(sumtypes==nC(find(nT==min(nT)))))); % indices of C that are from the conditions with min number of trials
            randi = randsample(1:max(nT), min(nT));
            cond_select = sort([maxi(randi);mini]);
            cond_reject = 1:D.ntrials;
            cond_reject(cond_select) = [];
            D = reject(D, cond_reject, 1)
            save(D)
        end

        % Average over trials per condition
        %==========================================================================
        clear S D;
        S.robust = 0;
        S.review = 0;
        %for i = 1:subspergrp
            fprintf('Averaging: subject %i\n',i);
            fname = ['a' fnames{i}];
            S.D = fullfile(pwd,fname)
            D = spm_eeg_average(S);
            %save(D);
        %end
%    end
%end

        % Convert to 3D images (2 space and 1 time dimension)
        %==========================================================================
        clear S D;
        S.images.fmt = 'channels';
        %S.images.elecs =
        %S.images.region_no =
        %S.images.freqs =
        %S.images.t_win =
        S.n = 32;
        S.interpolate_bad = 1;
        %for i = 1:subspergrp % [1:9 11:13]
        fprintf('Converting to images: subject %i\n',i);
        fname = ['ma' fnames{i}]
         %   load(fname);
         %   D.sensors.eeg.pnt = D.sensors.eeg.pnt(:,[2 1 3]);
         %   %D.sensors.eeg.elecpos = D.sensors.eeg.elecpos(:,[2 1 3]);
         %   D.fiducials.pnt = D.fiducials.pnt(:,[2 1 3]);
         %   D.fiducials.fid.pnt = D.fiducials.fid.pnt(:,[2 1 3]);
         %   D=meeg(D);
         %   S1 = [];
         %   S1.task = 'project3D';
         %   S1.modality = 'EEG';
         %   S1.updatehistory = 0;
         %   S1.D = D;
         %   D = spm_eeg_prep(S1);
         %   save(D);
        %end
    %end

     %   

            S.D = fullfile(pwd,fname);
            %S.mode = 'scalp x time';
            [D, S, Pout] = spm_eeg_convert2images(S) % spm8
            %images = spm_eeg_convert2images(S) % spm12
        %end



    %%
        % Smooth images
        %==========================================================================
        load('batch_smooth');
        matlabbatch{1,1}.spm.spatial.smooth.fwhm = [20 20 0]; % space space time
        matlabbatch{1,1}.spm.spatial.smooth.dtype = 0;
        matlabbatch{1,1}.spm.spatial.smooth.im = 1; % implicit mask
        %for is = start:length(fnames)
            fname = ['ma' fnames{i}];
            d = dir(fullfile(pwd,fname,'type*'));%spm12
            isub = [d(:).isdir]; %# returns logical vector%spm8
            nameFolds = {d(isub).name}';%spm8
            %nameFolds = {d.name}';%spm12
            nameFolds(ismember(nameFolds,{'.','..'})) = [];%spm8
            for nf = 1:length(nameFolds)
                imfile = dir(fullfile(pwd,fname,nameFolds{nf},'trial*.img'));%spm8
                %imfile = dir(fullfile(pwd,fname,nameFolds{nf}));%spm12
                imfile = imfile.name;
                matlabbatch{1,1}.spm.spatial.smooth.data = {fullfile(pwd,fname,nameFolds{nf},[imfile ',1'])};
                spm_jobman('initcfg')
                spm_jobman('run',matlabbatch);
            end
        %end
    end
end
    matlabmail
