function S=FT_grand_average(S)

dbstop if error

switch S.select.datatype
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

% unique indices of S.filelist to include in separate grand averages
col_ind = find(ismember(S.designmat(1,:),S.grand_avg.parts));
designtab = cell2table(S.designmat(2:end,col_ind)); % convert to table because unique with rows does not work on cell arrays!
[~,first_ind,file_ind]=unique(designtab,'rows','stable');
uni_ind = unique(file_ind);

% file loop
data_all = {}; % empty cell array for all subjects' tldata
S.gadata = {}; % empty cell array for grand average data
for ga = 1:length(uni_ind)
    
    % FIND THE FILES
    S.gafiles{ga} = S.filelist(file_ind==uni_ind(ga));
    
    % function to select events and frequencies
    S.files = S.gafiles{ga};
    S = select_eventsfreqs(S);
    data_all{ga} = S.data;
    
    method = {'across','within'}; if S.grand_avg.weighted; method = method{2}; else method = method{1};end
    data_all{ga}=horzcat(data_all{ga}{:});
    
    switch S.select.datatype
        case {'ERP'}
            % grand average using Fieldtrip
            cfg.channel        = data_all{ga}{1}.label(S.inclchan);
            cfg.latency        = 'all';
            cfg.keepindividual = 'no';
            cfg.normalizevar   = 'N-1';
            cfg.method         = method;
            cfg.parameter      = 'avg';
            S.gadata{ga} = ft_timelockgrandaverage_cab(cfg, data_all{ga});
        case {'Freq','TF'}
            cfg.keepindividual = 'no';
            cfg.foilim         = 'all'; %[fmin fmax] or 'all', to specify a subset of frequencies (default = 'all')
            cfg.toilim         = 'all'; % to specify a subset of latencies (default = 'all')
            cfg.channel        = data_all{ga}{1}.label(S.inclchan);
            cfg.parameter      = 'powspctrm';
            S.gadata{ga} = ft_freqgrandaverage_cab(cfg, data_all{ga});
    end
    
    gadata = S.gadata{ga};
    tabcell = table2cell(designtab(first_ind(ga),:));
    S.ganame{ga} = [strjoin(tabcell,'') 'grandavg'];
    switch S.select.datatype
        case {'TF','Freq'}
            fr_name = '_freq';
            for fr = 1:length(S.select.freqs)
                fr_name = [fr_name '_' num2str(S.select.freqs(fr))];
            end
            S.ganame{ga} = [S.ganame{ga} fr_name];
    end
    save(fullfile(S.path.file,S.ganame{ga}),'gadata');

end

