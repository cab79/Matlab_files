function h = CreateSequence(h)

disp('Creating sequence...');
h.Seq =struct;


%% create final sequences/blocks
if ~isfield(h.Seq,'signal')
    h.Seq.signal=1;
    h.Seq.condnum=1;
    h.Seq.blocks=1;
end