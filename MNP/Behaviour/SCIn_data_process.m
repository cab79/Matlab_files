function [S,D] = SCIn_data_process(S,D)
% requires strcuture D containing all the data

for d = 1:length(D)
    
    for op = 1:length(D(d).Output)

        % get settings if not in D already
        if ~isfield(D(d).Output(op),'Settings')
            A=load(fullfile('C:\Data\Matlab\Matlab_files\NTIP\SCIn\Sequences', D(d).Output(op).SeqName));
            D(d).Output(op).Settings = A.settings;
        end

        if S.accuracy.on
            % create presssignal 
            D(d).Processed(op).presssignal = nan(1,length(D(d).Sequence.signal));
            presssignal = [];
            for i = 1:length(S.accuracy.buttons)
                % translate actual buttons to their "signal" meaning
                presssignal(strcmp(D(d).Output(op).pressbutton,S.accuracy.buttons{i}))=S.accuracy.signal(i);
            end
            D(d).Processed(op).presssignal(D(d).Output(op).presstrial) = presssignal;

            % compare to actual signal
            D(d).Processed(op).correct = double(D(d).Sequence.signal==D(d).Processed(op).presssignal);
            D(d).Processed(op).correct(isnan(D(d).Processed(op).presssignal)) = nan;

            % split into conditions
            conds = unique(D(d).Sequence.condnum);
            for i = 1:length(conds)
                D(d).Processed(op).condcorrect{i} = D(d).Processed(op).correct(D(d).Sequence.condnum==conds(i)); 
                D(d).Processed(op).condcorrectfract(i) = nansum(D(d).Processed(op).condcorrect{i})/sum(~isnan(D(d).Processed(op).condcorrect{i}));
            end
        end
    end
    
end