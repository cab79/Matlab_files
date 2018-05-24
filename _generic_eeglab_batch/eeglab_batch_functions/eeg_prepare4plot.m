function S=eeg_prepare4plot(S)

dbstop if error
S.func = 'plot';
switch S.(S.func).select.datatype
    case 'ERP'
        S.path.file = S.path.erp;
    case 'TF'
        S.path.file = S.path.tf;
    case 'Freq'
        S.path.file = S.path.freq;
end

% GET FILE LIST
S = getfilelist(S);

% select channels
S=select_chans(S);

% unique indices of S.(S.func).filelist to include in separate grand averages
col_ind = find(ismember(S.(S.func).designmat(1,:),S.(S.func).parts));
designtab = cell2table(S.(S.func).designmat(2:end,col_ind)); % convert to table because unique with rows does not work on cell arrays!
[~,first_ind,file_ind]=unique(designtab,'rows','stable');
uni_ind = unique(file_ind);

data_all = {}; % empty cell array for all subjects' file data
S.(S.func).plotdat = {}; % empty cell array for data to plot
S.(S.func).plotmovdat = {}; % empty cell array for moving-averaged data

% file loop: for each unique plot
for pl = 1:length(uni_ind)
    
    % FIND THE FILES
    S.(S.func).plotfiles{pl} = S.(S.func).filelist(file_ind==uni_ind(pl));
    
    % function to select events and frequencies
    S.(S.func).files = S.(S.func).plotfiles{pl};
    S = select_eventsfreqs(S);
    data_all{pl} = S.(S.func).data;
    
    % average all the data
    data_all{pl}=horzcat(data_all{pl}{:}{:}); % first, make into a single struct
    switch S.(S.func).select.datatype
        case 'ERP'
            nd=ndims(data_all{pl}(1).avg);
            S.(S.func).plotdat{pl} = mean(cat(nd+1,data_all{pl}.avg),nd+1);
            %for d = 1:length(data_all{pl})
            %end
            
        case {'TF','Freq'}
            nd=ndims(data_all{pl}(1).powspctrm);
            S.(S.func).plotdat{pl} = mean(cat(nd+1,data_all{pl}.powspctrm),nd+1);
            if ~any(ismember(S.plot.parts,'freqs'));
                S.(S.func).plotdat{pl} = mean(S.(S.func).plotdat{pl},2);
            end
            % if freq, apply to moving average data
            if isfield(data_all{pl},'madata')
                nd=ndims(data_all{pl}(1).madata);
                S.(S.func).plotmovdat{pl} = mean(cat(nd+1,data_all{pl}.madata),nd+1);
                if ~any(ismember(S.plot.parts,'freqs'));
                    S.(S.func).plotmovdat{pl} = mean(S.(S.func).plotmovdat{pl},2);
                end
            end
    end
end

% Further data operations
if isfield(S.plot,'op')
    for oi = 1:length(S.plot.op)
        switch S.plot.op(oi).type
            case 'freq'
                opdat = S.(S.func).plotdat; 
                S.(S.func).plotdat = [];
                for c = 1:length(opdat)
                    for su = 1:size(S.plot.op(oi).grouping,1)
                        fni = 1:size(S.plot.op(oi).grouping,1)
                        fni(fi)=[];
                        S.(S.func).plotdat{c}(:,fi) = basenorm(opdat{c}(:,fi),mean(opdat{c}(:,fni),2),S.plot.op(oi).operation);
                    end
                    if isfield(S.(S.func),'plotmovdat')
                        opmovdat = S.(S.func).plotmovdat;
                        S.(S.func).plotmovdat = []; 
                        for su = 1:size(S.plot.op(oi).grouping,1)
                            fni = 1:size(S.plot.op(oi).grouping,1)
                            fni(fi)=[];
                            S.(S.func).plotmovdat{c}(:,fi,:) = basenorm(opmovdat{c}(:,fi,:),mean(opmovdat{c}(:,fni,:),2),S.plot.op(oi).operation);
                        end
                    end
                end


            case 'event'
                opdat = S.(S.func).plotdat; 
                S.(S.func).plotdat = [];
                for su = 1:size(S.plot.op(oi).grouping,1)
                    dat1 = opdat{S.plot.op(oi).grouping(su,1)};
                    dat2 = opdat{S.plot.op(oi).grouping(su,2)};
                    S.(S.func).plotdat{su} = basenorm(dat1,dat2,S.plot.op(oi).operation);
                end
                if isfield(S.(S.func),'plotmovdat')
                    opmovdat = S.(S.func).plotmovdat;
                    S.(S.func).plotmovdat = []; 
                    for su = 1:size(S.plot.op(oi).grouping,1)
                        dat1 = opmovdat{S.plot.op(oi).grouping(su,1)};
                        dat2 = repmat(mean(opmovdat{S.plot.op(oi).grouping(su,2)},3),1,1,size(dat1,3));
                        S.(S.func).plotmovdat{su}= basenorm(dat1,dat2,S.plot.op(oi).operation);
                    end
                end
        end
    end
end

function data=basenorm(data,meanVals,baselinetype)
if (strcmp(baselinetype, 'subtract'))
  data = data - meanVals;
elseif (strcmp(baselinetype, 'relative'))
  data = data ./ meanVals;
elseif (strcmp(baselinetype, 'relchange'))
  data = (data - meanVals) ./ meanVals;
elseif (strcmp(baselinetype, 'normchange')) || (strcmp(baselinetype, 'vssum'))
  data = (data - meanVals) ./ (data + meanVals);
elseif (strcmp(baselinetype, 'db'))
  data = 10*log10(data ./ meanVals);
else
  error('unsupported method for baseline normalization: %s', baselinetype);
end
