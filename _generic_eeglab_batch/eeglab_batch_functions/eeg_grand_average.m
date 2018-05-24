function S=FT_grand_average(S)

dbstop if error

S.func = 'ga';
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

% unique indices of S.filelist to include in separate grand averages
col_ind = find(ismember(S.(S.func).designmat(1,:),S.(S.func).grand_avg.parts));
designtab = cell2table(S.(S.func).designmat(2:end,col_ind)); % convert to table because unique with rows does not work on cell arrays!
[~,first_ind,file_ind]=unique(designtab,'rows','stable');
uni_ind = unique(file_ind);

% file loop
data_all = {}; % empty cell array for all subjects' tldata
S.(S.func).gadata = {}; % empty cell array for grand average data
for ga = 1:length(uni_ind)
    
    % FIND THE FILES
    S.(S.func).gafiles{ga} = S.(S.func).filelist(file_ind==uni_ind(ga));
    
    % function to select events and frequencies
    S.(S.func).files = S.(S.func).gafiles{ga};
    S = select_eventsfreqs(S);
    data_all{ga} = S.(S.func).data;
    
    method = {'across','within'}; if S.(S.func).grand_avg.weighted; method = method{2}; else method = method{1};end
    
    % get index of files each condition came from
    Nconds = cellfun(@size,data_all{ga},'UniformOutput',0);
    fileidx=[];
    for n = 1:length(Nconds)
        fileidx(n,1:max(Nconds{n})) = n;
    end
    fileidx = reshape(fileidx',1,[]);
    data_all{ga}=horzcat(data_all{ga}{:});
    
    % remove empty cells
    noemp = find(~cellfun(@isempty,data_all{ga}));
    data_all{ga} = data_all{ga}(noemp);
    S.(S.func).fileidx = fileidx(noemp);
    
    % savename
    tabcell = table2cell(designtab(first_ind(ga),:));
    S.(S.func).ganame{ga} = [strjoin(tabcell,'') 'grandavg'];
    switch S.(S.func).select.datatype
        case {'TF','Freq'}
            fr_name = '_freq';
            for fr = 1:length(S.(S.func).select.freqs)
                fr_name = [fr_name '_' num2str(S.(S.func).select.freqs(fr))];
            end
            S.(S.func).ganame{ga} = [S.(S.func).ganame{ga} fr_name];
    end
    
    % multivariate outliers
    if S.ga.grand_avg.outliers
        S=MultiOutliers(S,data_all{ga});
    end
    outdata = S.(S.func).multout;
    outlist = S.(S.func).multoutlist;
    save(fullfile(S.path.file,'Outliers.mat'),'outdata','outlist');
    
    switch S.(S.func).select.datatype
        case {'ERP'}
            % grand average using Fieldtrip
            cfg.channel        = data_all{ga}{1}.label(S.(S.func).inclchan);
            cfg.latency        = 'all';
            cfg.keepindividual = 'no';
            cfg.normalizevar   = 'N-1';
            cfg.method         = method;
            cfg.parameter      = 'avg';
            S.(S.func).gadata{ga} = ft_timelockgrandaverage_cab(cfg, data_all{ga});
        case {'Freq','TF'}
            cfg.keepindividual = 'no';
            cfg.foilim         = 'all'; %[fmin fmax] or 'all', to specify a subset of frequencies (default = 'all')
            cfg.toilim         = 'all'; % to specify a subset of latencies (default = 'all')
            cfg.channel        = data_all{ga}{1}.label(S.(S.func).inclchan);
            cfg.parameter      = 'powspctrm';
            S.(S.func).gadata{ga} = ft_freqgrandaverage_cab(cfg, data_all{ga});
    end
    
    gadata = S.(S.func).gadata{ga};
    save(fullfile(S.path.file,S.(S.func).ganame{ga}),'gadata');
    

end

