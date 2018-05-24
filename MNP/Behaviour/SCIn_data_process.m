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
            D(d).Processed(op).presssignal = nan(1,size(D(d).Sequence.signal,2));
            presssignal = [];
            for i = 1:length(S.accuracy.buttons)
                % translate actual buttons to their "signal" meaning
                presssignal(strcmp(D(d).Output(op).pressbutton,S.accuracy.buttons{i}))=S.accuracy.signal(i);
            end
            D(d).Processed(op).presssignal(D(d).Output(op).presstrial) = presssignal;
            
            % what are the correct responses?
            D(d).Processed(op).correct_resp = D(d).Sequence.signal(S.signal.target,:);
            D(d).Processed(op).correct_resp_new = D(d).Processed(op).correct_resp;
            for i = 1:length(S.accuracy.target_resp{1})
                D(d).Processed(op).correct_resp_new(D(d).Processed(op).correct_resp==S.accuracy.target_resp{1}(i)) = S.accuracy.target_resp{2}(i);
            end
            D(d).Processed(op).correct_resp = D(d).Processed(op).correct_resp_new;
            
            % compare actual responses to correct responses - which were
            % correct?
            D(d).Processed(op).correct = double(D(d).Processed(op).correct_resp==D(d).Processed(op).presssignal);
            D(d).Processed(op).correct(isnan(D(d).Processed(op).presssignal)) = nan;
            
            % moving average percent correct over time
            ma=S.movingavg;
            for i = 1:length(D(d).Processed(op).correct)
                D(d).Processed(op).macorrect(i) = 100*sum(D(d).Processed(op).correct(max(1,i-ma+1):i))/length(max(1,i-ma+1):i);
            end
            
            % create condnum cell which either contains all trials or a
            % subset of trials according to S.trialmax
                    % This allows testing of how many trials are needed to
                    % produce robust condition differences.
            condnum={};
            conds = unique(D(d).Sequence.condnum);
            if ~isempty(S.trialmax)
                for tm = 1:length(S.trialmax)
                    trialmax = S.trialmax{tm};
                    % reduce the number of trials per condition-pair to S.trialmax
                    for i = 1:length(conds)/2
                        ind = (i-1)*2+1:(i-1)*2+2; % condition pair indices, e.g. 1 2, 3 4, 5 6
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

            % for each version of condnum (diff number of trials in each),
            % get perc correct for each condition and block
            for cn = 1:length(condnum)
                
                % split into conditions
                condsuni = unique(condnum{cn});
                condsuni = condsuni(~isnan(condsuni));
                D(d).Processed(op).condcorrectfract{cn} = nan(1,length(conds));
                for i = 1:length(condsuni)
                    D(d).Processed(op).condcorrect{cn}{condsuni(i)} = D(d).Processed(op).correct(condnum{cn}==condsuni(i) & ismember(D(d).Sequence.signal(S.signal.target,:),S.accuracy.cond_stimtype)); 
                    D(d).Processed(op).numtrials{cn}{condsuni(i)} = sum(~isnan(D(d).Processed(op).condcorrect{cn}{condsuni(i)}));
                    D(d).Processed(op).condcorrectfract{cn}(condsuni(i)) = nansum(D(d).Processed(op).condcorrect{cn}{condsuni(i)})/D(d).Processed(op).numtrials{cn}{condsuni(i)};
                end

                % split into blocks
                blocks = unique(D(d).Sequence.blocks);
                for i = 1:length(condsuni)
                    for b = 1:length(blocks)
                        D(d).Processed(op).blockcondcorrect{cn}{condsuni(i)}{b} = D(d).Processed(op).correct(condnum{cn}==condsuni(i) & D(d).Sequence.blocks==blocks(b) & ismember(D(d).Sequence.signal(S.signal.target,:),S.accuracy.cond_stimtype)); 
                        D(d).Processed(op).blocknumtrials{cn}{condsuni(i)}{b} = sum(~isnan(D(d).Processed(op).blockcondcorrect{cn}{condsuni(i)}{b}));
                        D(d).Processed(op).blockcondcorrectfract{cn}{condsuni(i)}(b) = nansum(D(d).Processed(op).blockcondcorrect{cn}{condsuni(i)}{b})/D(d).Processed(op).blocknumtrials{cn}{condsuni(i)}{b};
                        for ii = 1:length(D(d).Processed(op).blockcondcorrect{cn}{condsuni(i)}{b})
                            D(d).Processed(op).blockcondcorrectmovavg{cn}{condsuni(i)}{b}(ii) = 100*sum(D(d).Processed(op).blockcondcorrect{cn}{condsuni(i)}{b}(max(1,ii-ma+1):ii))/length(max(1,ii-ma+1):ii);;
                        end
                    end
                end
                
                
                % split into blocks before conditions, to get moving average of condition ratio of
                % percent correct
                blocks = unique(D(d).Sequence.blocks);
                for b = 1:length(blocks)
                    trialidx=find(D(d).Sequence.blocks==blocks(b) & ismember(D(d).Sequence.signal(S.signal.target,:),S.accuracy.cond_stimtype));
                    D(d).Processed(op).blockcorrect{cn}{b} = D(d).Processed(op).correct(trialidx);
                    cond = condnum{cn}(trialidx);
                    % for each trial, get indices of last ma trials
                    for ii = 1:length(trialidx)
                        idx=max(1,ii-ma+1):ii;
                        correcttrials=D(d).Processed(op).blockcorrect{cn}{b}(idx);
                        condidx=cond(idx);
                        ucond = unique(condidx);
                        correct=[];
                        Ntrials=[];
                        for ci =1:length(ucond)
                            % proportion correct, from 0 to 1: correctN / totalN
                            Ntrials(ci)=length(idx(condidx==ucond(ci)));
                            correct(ci) = nansum(D(d).Processed(op).blockcorrect{cn}{b}(idx(condidx==ucond(ci))))/Ntrials(ci);
                        end
                        if isempty(correct) || length(correct)<2 || any(Ntrials<2); 
                            D(d).Processed(op).blockcorrectmovavg{cn}{b}(ii) = NaN;
                        else
                            % difference in proportion correct ranging from -1 to 1
                            D(d).Processed(op).blockcorrectmovavg{cn}{b}(ii) = correct(1)-correct(2);
                        end
                    end
                end
                
                % split into stim intensity per condition
                stims = unique(D(d).Sequence.signal(S.signal.target,:));
                for i = 1:length(condsuni)
                    for s = 1:length(stims)
                        D(d).Processed(op).stimcondcorrect{cn}{condsuni(i)}{s} = D(d).Processed(op).correct(condnum{cn}==condsuni(i) & D(d).Sequence.signal(S.signal.target,:)==stims(s)); 
                        D(d).Processed(op).stimnumtrials{cn}{condsuni(i)}{s} = sum(~isnan(D(d).Processed(op).stimcondcorrect{cn}{condsuni(i)}{s}));
                        D(d).Processed(op).stimcondcorrectfract{cn}{condsuni(i)}(s) = nansum(D(d).Processed(op).stimcondcorrect{cn}{condsuni(i)}{s})/D(d).Processed(op).stimnumtrials{cn}{condsuni(i)}{s};
                    end
                end
                
                % split into stim intensity per cue
                if S.signal.cue
                    cues = unique(D(d).Sequence.signal(S.signal.cue,:));
                    for i = 1:length(cues)
                        for s = 1:length(stims)
                            D(d).Processed(op).stimcuecorrect{cn}{cues(i)}{s} = D(d).Processed(op).correct(D(d).Sequence.signal(S.signal.cue,:)==cues(i) & D(d).Sequence.signal(S.signal.target,:)==stims(s)); 
                            D(d).Processed(op).stimcuenumtrials{cn}{cues(i)}{s} = sum(~isnan(D(d).Processed(op).stimcuecorrect{cn}{cues(i)}{s}));
                            D(d).Processed(op).stimcuecorrectfract{cn}{cues(i)}(s) = nansum(D(d).Processed(op).stimcuecorrect{cn}{cues(i)}{s})/D(d).Processed(op).stimcuenumtrials{cn}{cues(i)}{s};
                        end
                    end
                end
            end
        end
    end
    
end