cfg = [];
cfg.layout = 'acticap-64ch-standard2.mat';
cfg.ylim = [-15 15]
for d = 1:length(tldata)
    cfg.dataname{d} = ['cond' num2str(d)];
end

f=figure('units','normalized','outerposition',[0 0 1 1]);
ft_multiplotER_cab(cfg, tldata);
title('title')

% good, bad?
qual{i,1} = fname;
qual{i,2} = questdlg('Data quality','','Good','So-so','Bad','So-so')
close(f)