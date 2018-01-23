function h = SimpleSequence(h)

disp('Creating sequence...');
h.Seq =struct;


%% create final sequences/blocks
if ~isfield(h.Seq,'signal')
    h.Seq.signal=ones(1,h.Settings.ntrials);
    h.Seq.condnum=ones(1,h.Settings.ntrials);
    h.Seq.blocks=ones(1,h.Settings.ntrials);
end