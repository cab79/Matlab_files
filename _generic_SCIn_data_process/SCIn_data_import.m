function [S,D] = SCIn_data_import(S)
% takes outputs from SCIn and puts them in a structure 'D', with field name
% being the S.loadprefixes

% get filename using only first input type (filename prefix)
S.load.prefix = S.load.prefixes(1);

% GET FILE LIST
S.path.file = S.path.raw;
S = getfilelist(S);

% change to the input directory
eval(sprintf('%s', ['cd(''' S.path.file ''')']));

% report if there are no such files
if isempty(S.filelist)
    error('No files found to import!\n');
end

% run though all files in a loop and get data
for s = 1:length(S.select.subjects)
    
    % create output D
    D(s).subname = S.select.subjects{s};
        
    % FIND THE FILES FOR THIS SUBJECT
    subfiles = S.filelist(find(not(cellfun('isempty', strfind(S.filelist,D(s).subname)))));
    dirnames = S.dirlist(find(not(cellfun('isempty', strfind(S.filelist,D(s).subname)))));

    % CYCLE THROUGH EACH FILE FOR THIS SUBJECT
    for f = 1:length(subfiles)
        filename = subfiles{f};
        dirname = dirnames{f};
        S.load.prefix = S.load.prefixes{1};
        fprintf('\nProcessing %s.\n\n', filename);

        % assume only one field in loaded file
        temp = load(fullfile(dirname,filename));
        fename = fieldnames(temp(f));
        if isstruct(temp.(fename{1}))
            outstruct = temp.(fename{1});
        else
            outstruct.Output = temp.(fename{1});
        end
        D(s).(S.load.prefix)(f) = outstruct;

        % get other filetypes with same name
        for fn = 1:length(S.load.prefixes)-1
            S.load.prefix = S.load.prefixes{fn+1};
            filename = strrep(filename,S.load.prefixes{1},S.load.prefix);

            % assume only one field in loaded file
            temp = load(fullfile(dirname,filename));
            fename = fieldnames(temp(f));
            D(s).(S.load.prefix)(f) = temp.(fename{1});
        end
    end
end