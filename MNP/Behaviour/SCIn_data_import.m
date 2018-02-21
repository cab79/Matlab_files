function [S,D] = SCIn_data_import(S)
% takes outputs from SCIn and puts them in a structure 'D', with field name
% being the S.loadprefixes

% get filename using only first input type (filename prefix)
S.loadprefix = S.loadprefixes{1};

% GET FILE LIST
S.filepath = S.rawpath;
S = getfilelist(S);

% change to the input directory
eval(sprintf('%s', ['cd(''' S.rawpath ''')']));

% report if there are no such files
if isempty(S.filelist)
    error('No files found to import!\n');
end

% run though all files in a loop and get data
for s = 1:length(S.subjects)
    
    % create output D
    D(s).subname = S.subjects{s};
        
    % FIND THE FILES FOR THIS SUBJECT
    subfiles = S.filelist(find(not(cellfun('isempty', strfind(S.filelist,D(s).subname)))));

    % CYCLE THROUGH EACH FILE FOR THIS SUBJECT
    for f = 1:length(subfiles)
        filename = subfiles{f};
        S.loadprefix = S.loadprefixes{1};
        fprintf('\nProcessing %s.\n\n', filename);

        % assume only one field in loaded file
        temp = load(filename);
        fename = fieldnames(temp(f));
        D(s).(S.loadprefix)(f) = temp.(fename{1});

        % get other filetypes with same name
        for fn = 1:length(S.loadprefixes)-1
            S.loadprefix = S.loadprefixes{fn+1};
            filename = strrep(filename,S.loadprefixes{1},S.loadprefix);

            % assume only one field in loaded file
            temp = load(filename);
            fename = fieldnames(temp(f));
            D(s).(S.loadprefix)(f) = temp.(fename{1});
        end
    end
end