clear all

restoredefaultpath
run('C:\Data\Matlab\Matlab_files\CORE\CORE_addpaths')
%% SPECIFY DATA
filepath = 'C:\Data\CORE\EEG\ana\spm\SPMdata\sensorimages'; 
datfile = ['C:\Data\CORE\participants\Participant_data.xlsx']; % .xlsx file to group participants; contains columns named 'Subject', 'Group', and any covariates of interest

% prefix, middle part, or suffix of files to load (or leave empty) to select a subset of folders
%fpref = 't-200_899_b-200_0_mspm12_fnums';
%fpref = 't-200_899_b-200_0_mspm12';
fpref = 't-200_299_b-200_0_mspm12_fnums_flip';
fmid = '';
%fsuff = '_4_cleaned_tm';
%fsuff = '_2_merged_cleaned';
fsuff = '_4_merged_cleaned'; % fnum

% conditions to flip:
fcond = {[1:4],[5:8]}; % part4 fnum

%% RUN
fname=[fpref '*' fmid  '*' fsuff];
fname=strrep(fname,'**','*');
files = dir(fullfile(filepath,fname));

% get side affected data from participant file
subhead = 'Subject';
grphead = 'Group';
inchead = 'Include';
CRPSsidehead = 'CRPSlr';
include_codes = [1];
[~,~,pdata] = xlsread(datfile);
grp_col = find(strcmp(pdata(1,:),grphead));
sub_col = find(strcmp(pdata(1,:),subhead));
inc_col = find(strcmp(pdata(1,:),inchead));
side_col = find(strcmp(pdata(1,:),CRPSsidehead));
inc_idx = cellfun(@(x) ismember(x,include_codes), pdata(2:end,inc_col), 'UniformOutput', 0);
inc_idx = find(cell2mat(inc_idx));
subjects = pdata(2:end,sub_col);
CRPSsides = pdata(2:end,side_col);

% add random sides to the healthy subjects (who have no CRPS side)
sidenan = cell2mat(cellfun(@(x) any(x), cellfun(@(x) isnan(x), CRPSsides, 'UniformOutput', 0), 'UniformOutput', 0));
sidepool = CRPSsides(~sidenan);
sn_idx=find(sidenan);
for s=1:length(sn_idx)
    CRPSsides(sn_idx(s)) = sidepool(randi([1 numel(sidepool)]));
end

allfiles = {files(:).name};
for f = 1:length(subjects)
    
    try
        fname = allfiles{~cellfun(@isempty,strfind(allfiles,subjects{f}))};
    catch
        continue
    end
 
    % swap the images
    for nf = 1:length(fcond{1})
        input1 = fullfile(filepath,fname,['scondition_' num2str(fcond{1}(nf)) '.nii']);
        input2 = fullfile(filepath,fname,['scondition_' num2str(fcond{2}(nf)) '.nii']);
        
        if strcmp(CRPSsides{f},'R')
            output2 = fullfile(filepath,fname,['scondition_' num2str(fcond{1}(nf)) '_swapped.nii']);
            output1 = fullfile(filepath,fname,['scondition_' num2str(fcond{2}(nf)) '_swapped.nii']);
        else
            output1 = fullfile(filepath,fname,['scondition_' num2str(fcond{1}(nf)) '_swapped.nii']);
            output2 = fullfile(filepath,fname,['scondition_' num2str(fcond{2}(nf)) '_swapped.nii']);
        end
        
        copyfile(input1,output1);
        copyfile(input2,output2);
        delete(input1);
        delete(input2);
    end
end
