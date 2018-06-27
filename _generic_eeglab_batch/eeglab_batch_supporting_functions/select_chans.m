function S=select_chans(S);

if ~isfield(S.(S.func),'chanlocs')
    load(fullfile(S.path.main,'chanlocs.mat'));
    S.(S.func).chanlocs = chanlocs;
end
if isempty(S.(S.func).select.chans{1})
    S.(S.func).inclchan = 1:length(S.(S.func).chanlocs);
else
    S.(S.func).inclchan = S.(S.func).select.chans{1};
end
if ~isempty(S.(S.func).select.chans{2})
    S.(S.func).inclchan(S.(S.func).select.chans{2}) = [];
end