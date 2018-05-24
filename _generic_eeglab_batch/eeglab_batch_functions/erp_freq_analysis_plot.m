function S=erp_freq_analysis_plot(S)

%% get subject file info and load

% GET FILE LIST
S.filepath = S.freqpath;
S.loadext='Freq.mat';
S = getfilelist(S);

% select session, block, frequency
S.fn = dsearchn(S.freqsrange',S.freqselect);

% select channels
if ~exist('chanlocs','var')
    load(fullfile(S.setpath,'chanlocs.mat'));
    S.chanlocs = chanlocs;
end
S.inclchan = 1:length(S.chanlocs);
S.inclchan(S.exclchan) = [];

% load subject data
for s = 1:length(S.subjects)
    %if participant_eventorder(s)==-1
    %    event_order = fliplr(event_reorder);
    %else
    %    event_order = event_reorder;
    %end
    
    % FIND THE FILES FOR THIS SUBJECT
    subfiles = S.filelist(find(not(cellfun('isempty', strfind(S.filelist,S.subjects{s})))));
    % FIND THE FILES FOR THIS BLOCK
    blockfiles = subfiles(find(not(cellfun('isempty', strfind(subfiles,S.blocks{S.blockselect})))));
    
    % get the data
    for c = 1:length(S.conds)
        % get filename
        file = blockfiles(find(not(cellfun('isempty', strfind(blockfiles,S.conds{c})))));
        if length(file)~=1
            error('file name not uniquely specified')
        end
        % load data
        load(fullfile(S.filepath,file{:}));
        % average data over selected conditions
        if isstruct(fdata)
            dat = fdata.powspctrm; %trials,chans,freqs
            madat = fdata.madata;
            act_freq = fdata.freq;
        elseif iscell(fdata)
            clear datc madatc
            for e = 1:length(S.eventselect)
                datc{e} = fdata{S.eventselect(e)}.powspctrm;
                madatc{e} = fdata{S.eventselect(e)}.madata;
            end
            dat = cat(1,datc{:});
            madat = cat(4,madatc{:});
            act_freq = fdata{1}.freq;
        end
        S.meandat{c,1} = squeeze(mean(dat,1));
        S.movmeandat{c,1} = squeeze(mean(madat,4));
    end

    if ~all(ismember(S.freqsrange,act_freq))
        error('output frequencies do not match the settings specified')
    end


    S.plotdat=S.meandat;
    S.plotmovdat=S.movmeandat;
    


end

