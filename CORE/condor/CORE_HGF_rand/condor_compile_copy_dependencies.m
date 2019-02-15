%% SETTTINGS
clear all
restoredefaultpath
run('C:\Data\Matlab\Matlab_files\CORE\CORE_addpaths')
main_path = 'C:\Data\Matlab\Matlab_files\CORE\condor';

% function name to find dependencies
fnames = {
    'CORE_condor_fit_step1';
    'CORE_condor_fit_step2';
    'CORE_condor_fit_job';
    'CORE_condor_monitor_build';
    };

% folder names to include all files, regardless of whether they are listed
% in the function as dependencies, needed if within eval functions
fsname = {
    'C:\Data\Matlab\Matlab_files\_generic_HGF'
    'C:\Data\Matlab\HGF\HGFv5.0'
    'C:\Data\Matlab\Matlab_files\CORE\behaviour\GBM_configs'
    };

% name of folder to create to copy the files to
sname = 'CORE_HGF_rand';

%% RUN

% create save dir
pth = fullfile(main_path,sname);
if ~exist(pth,'dir')
    mkdir(pth)
end
cd(pth)

% run
fList = {};
for fn = 1:length(fnames)
    [fList_new,pList] = matlab.codetools.requiredFilesAndProducts([fnames{fn} '.m']);
    fList = [fList;fList_new'];
end

% get further dependencies and add to the end of the list
% f = 0;
% len=length(fList);
% while f<len
%     f = f+1;
%     len=length(fList);
%     disp(['adding dependency ' num2str(f) '/' num2str(len)])
%     [fList_sub,~] = matlab.codetools.requiredFilesAndProducts(fList{f});
%     
%     % identify unique additions
%     new_ones = ~ismember(fList_sub, fList);
%     
%     % add them on
%     fList = [fList;fList_sub(new_ones)'];
% end

for f = 1:length(fList)
    try
        disp(['copying file ' num2str(f) '/' num2str(length(fList))])
        fFullPath = fullfile(pth,strrep(fList{f},'C:\Data\Matlab','dependencies'));
        fDir = fileparts(fFullPath);
        if ~exist(fDir,'dir')
            mkdir(fDir)
        end
        copyfile(fList{f},fFullPath);
    catch
        disp(['could not copy ' fList{f}])
    end
end

% further folders (not obvious dependencies)
for f = 1:length(fsname)
    fFullPath = fullfile(pth,strrep(fsname{f},'C:\Data\Matlab','dependencies'));
    copyfile(fsname{f},fFullPath)
end

% create zip and delete unzipped folder
zip('dependencies.zip',fullfile(pth,'dependencies'))
eval(['rmdir ' fullfile(pth,'dependencies') ' s']);