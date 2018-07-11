function [S,D] = SCIn_seq_import(S)

S.path.file = S.path.seq;
eval(sprintf('%s', ['cd(''' S.path.file ''')']));
load(S.select.blocks{1});

% create output D
D(1).subname = 'Seq';
D(1).Sequence = seq;
D(1).Output.Settings.buttonopt = {'LeftArrow', 'RightArrow'};
