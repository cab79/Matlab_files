function S = select_eventsfreqs(S)

S.(S.func).data={};
for f = 1:length(S.(S.func).files)
    
    disp(['Select data: file ' num2str(f) '/' num2str(length(S.(S.func).files))]);

    % load data for each file
    load(fullfile(S.path.file,S.(S.func).files{f}));

    switch S.(S.func).select.datatype
        case 'ERP'
            % select events
            if isstruct(tldata)
                S.(S.func).data{f} = {tldata};
            elseif iscell(tldata)
                S.(S.func).data{f} = tldata(S.(S.func).select.events);
            end
        case {'TF','Freq'}
            % select events
            if isstruct(fdata)
                S.(S.func).data{f} = {fdata};
            elseif iscell(fdata)
                S.(S.func).data{f} = fdata(S.(S.func).select.events);
            end
            % create subject averages for specified freq
            for ev = 1:length(S.(S.func).data{f})
                %S.(S.func).fn = dsearchn(S.(S.func).data{f}{ev}.freq',S.(S.func).select.freqs);
                cfg.variance      = 'no';%'yes' or 'no', estimate standard error in the standard way (default = 'no')
                cfg.jackknife     = 'no';%'yes' or 'no', estimate standard error by means of the jack-knife (default = 'no')
                cfg.keeptrials    = 'no';%'yes' or 'no', estimate single trial power (useful for fourier data) (default = 'no')
                cfg.channel       = 'all';%%Nx1 cell-array with selection of channels (default = 'all'),see FT_CHANNELSELECTION for details
                cfg.trials        = 'all';%'all' or a selection given as a 1xN vector (default = 'all')
                cfg.frequency     = [min(S.(S.func).select.freqs),max(S.(S.func).select.freqs)]; %[fmin fmax] or 'all', to specify a subset of frequencies (default = 'all')
                cfg.latency       = 'all';%[tmin tmax] or 'all', to specify a subset of latencies (default = 'all')
                temp = S.(S.func).data{f}{ev};
                S.(S.func).data{f}{ev} = ft_freqdescriptives(cfg, temp);
                if isfield(temp,'madata')
                    S.(S.func).data{f}{ev}.madata = temp.madata;
                end
            end
            
    end

end