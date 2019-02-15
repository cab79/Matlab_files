clear all
% badly written script - prone to errors due to try catch statments - REVISE!

restoredefaultpath
run('C:\Data\Matlab\Matlab_files\CORE\CORE_addpaths')
%% SPECIFY DATA
filepath = 'C:\Data\CORE\EEG\ana\spm\SPMdata\sensorimages'; 

% prefix, middle part, or suffix of files to load (or leave empty) to select a subset of folders
%fpref = 't-200_899_b-200_0_mspm12_CPavg';
fpref = 't-200_899_b-200_0_mspm12';
%fpref = 't-200_899_b-200_0_mspm12_fnums';
%fpref = 't-200_299_b-200_0_mspm12_fnums';
fmid = '';
%fsuff = '_4_cleaned_tm';
fsuff = '_2_cleaned_tm';
%fsuff = '_4_merged_cleaned'; % fnum

fcond = {
    [1 5]; % CP10, Odd, DC1
    [2 6]; % CP10, Odd, DC3
    [3 7]; % CP10, Stan, DC1
    [4 8]; % CP10, Stan, DC3
    [9 13]; % CP30, Odd, DC1
    [10 14]; % CP30, Odd, DC3
    [11 15]; % CP30, Stan, DC1
    [12 16]; % CP30, Stan, DC3
    [17 21]; % CP50, Odd, DC1
    [18 22]; % CP50, Odd, DC3
    [19 23]; % CP50, Stan, DC1
    [20 24]; % CP50, Stan, DC3
};


%% RUN
fname=[fpref '*' fmid  '*' fsuff];
fname=strrep(fname,'**','*');
files = dir(fullfile(filepath,fname));

for f = 6:44%:length(files)
    fname = files(f).name;
 
    for nf = 9:length(fcond)
        for nc = 1:length(fcond{nf,:})
            try
                fs = fullfile(filepath,files(f).name,['scondition_' num2str(fcond{nf}(nc)) '.nii']);
                if ~exist(fs,'file')
                    error('error')
                end
            catch
                fs = fullfile(filepath,files(f).name,['scondition_' num2str(fcond{nf}(nc)) '_flip.nii']);
            end
            inputnames{nc} = fs;
        end
        
        try
            sname=fullfile(filepath,files(f).name,['scondition_' num2str(nf) '.nii']);
            spm_imcalc_ui(inputnames,sname,'mean(X)',{1}); 
        catch
            try
                sname=fullfile(filepath,files(f).name,['scondition_' num2str(nf) '.nii']);
                copyfile(inputnames{1},sname); 
            catch
                sname=fullfile(filepath,files(f).name,['scondition_' num2str(nf) '.nii']);
                copyfile(inputnames{2},sname); 
            end
        end
        
    end
    
    % cleanup
    cl = [fcond{:}];
    cl=cl(:);
    for cf = 1:max(cl)
        try
            fs = fullfile(filepath,files(f).name,['scondition_' num2str(cf) '_flip.nii']);
            delete(fs)
        end
    end
    for cf = nf+1:max(cl)
        try
            fs = fullfile(filepath,files(f).name,['scondition_' num2str(cf) '.nii']);
            delete(fs)
        catch
            try
                fs = fullfile(filepath,files(f).name,['scondition_' num2str(cf) '_flip.nii']);
                delete(fs)
            end
        end
    end
    
end
