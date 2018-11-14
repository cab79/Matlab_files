clear all

dbstop if error
restoredefaultpath
run('C:\Data\Matlab\Matlab_files\CORE\CORE_addpaths')
%% SPECIFY DATA
filepath = 'C:\Data\CORE\EEG\ana\spm\SPMdata\sourceimages_GS'; 
datfile = ['C:\Data\CORE\participants\Participant_data.xlsx']; % .xlsx file to group participants; contains columns named 'Subject', 'Group', and any covariates of interest


% prefix, middle part, or suffix of files to load (or leave empty) to select a subset of folders
%fpref = 't-200_899_b-200_0_mspm12';
fpref = 'mspm12_fnums';
fmid = '';
%fsuff = '_4_cleaned_tm';
%fsuff = '_2_merged_cleaned';
fsuff = '_4_merged_cleaned'; % fnum

% timewins to flip
tcond = [1 2];

% conditions to flip
fcond = {[1:4],[5:8]}; % part4 fnum

%% RUN
fbase=[fpref '*' fmid  '*' fsuff];
fbase=strrep(fbase,'**','*');

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

for nt = 1:length(tcond)

    for nf = 1:length(fcond{1})
        files = dir(fullfile(filepath,[fbase '_' num2str(tcond(nt)) '*' num2str(fcond{1}(nf)) '.nii']));
        allfiles = {files(:).name};
        for f = 1:length(files)
            subind=[];
            for i = 1:length(subjects)
                if ~isempty(strfind(allfiles{f},subjects{i}))
                    subind=i;
                end
            end
            
            fname = fullfile(filepath,files(f).name); 
            [pth, nme, ext] = fileparts(fname);
            
            input1 = fname;
            input2 = strrep(fname,[num2str(fcond{1}(nf)) '.nii'],[num2str(fcond{2}(nf)) '_flip.nii']);

            if strcmp(CRPSsides{subind},'R')
                output2 = strrep(fname,[num2str(fcond{1}(nf)) '.nii'],[num2str(fcond{1}(nf)) '_swapped.nii']);
                output1 = strrep(fname,[num2str(fcond{1}(nf)) '.nii'],[num2str(fcond{2}(nf)) '_flip_swapped.nii']);
            else
                output1 = strrep(fname,[num2str(fcond{1}(nf)) '.nii'],[num2str(fcond{1}(nf)) '_swapped.nii']);
                output2 = strrep(fname,[num2str(fcond{1}(nf)) '.nii'],[num2str(fcond{2}(nf)) '_flip_swapped.nii']);
            end

            copyfile(input1,output1);
            copyfile(input2,output2);
            delete(input1);
            delete(input2);
        end

    end
end
