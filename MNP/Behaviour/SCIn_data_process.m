function [S,D] = SCIn_data_process(S,D)
% requires strcuture D containing all the data

for d = 1:length(D)
    
    for op = 1:length(D(d).Output)

        % get settings if not in D already
        if ~isfield(D(d).Output(op),'Settings')
            A=load(fullfile('C:\Data\Matlab\Matlab_files\NTIP\SCIn\Sequences', D(d).Output(op).SeqName));
            D(d).Output(op).Settings = A.settings;
        end
        
        S.accuracy.buttons = D(d).Output.Settings.buttonopt;

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
            
            condnum={};
            conds = unique(D(d).Sequence.condnum);
            if ~isempty(S.trialmax)
                for tm = 1:length(S.trialmax)
                    trialmax = S.trialmax{tm};
                    % reduce the number of trials per condition-pair to S.trialmax
                    
                    for i = 1:length(conds)/2
                        ind = (i-1)*2+1:(i-1)*2+2; % condition pair indices
                        trialind = find(ismember(D(d).Sequence.condnum,conds(ind)));
                        trialind_new = trialind(1:min(trialmax,length(trialind)));
                        condnum_orig{tm} = D(d).Sequence.condnum;
                        condnum{tm}(trialind) = nan;
                        condnum{tm}(trialind_new) = condnum_orig{tm}(trialind_new);
                    end
                end
            else
                condnum{1} = D(d).Sequence.condnum;
            end

            % for each version of condnum (diff number of trials in each)
            for cn = 1:length(condnum)
                % split into conditions
                condsuni = unique(condnum{cn});
                condsuni = condsuni(~isnan(condsuni));
                
                D(d).Processed(op).condcorrectfract{cn} = nan(1,length(conds));
                for i = 1:length(condsuni)
                    D(d).Processed(op).condcorrect{cn}{condsuni(i)} = D(d).Processed(op).correct(condnum{cn}==condsuni(i)); 
                    D(d).Processed(op).condcorrectfract{cn}(condsuni(i)) = nansum(D(d).Processed(op).condcorrect{cn}{condsuni(i)})/sum(~isnan(D(d).Processed(op).condcorrect{cn}{condsuni(i)}));
                end

                % split into blocks
                blocks = unique(D(d).Sequence.blocks);
                for i = 1:length(condsuni)
                    for b = 1:length(blocks)
                        D(d).Processed(op).blockcondcorrect{cn}{condsuni(i)}{b} = D(d).Processed(op).correct(condnum{cn}==condsuni(i) & D(d).Sequence.blocks==blocks(b)); 
                        D(d).Processed(op).blockcondcorrectfract{cn}{condsuni(i)}{b} = nansum(D(d).Processed(op).blockcondcorrect{cn}{condsuni(i)}{b})/sum(~isnan(D(d).Processed(op).blockcondcorrect{cn}{condsuni(i)}{b}));
                    end
                end
            end
        end
    end
    
end