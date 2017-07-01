%% plot FT data using FT functions
cfg = [];
%cfg.baseline     = [-0.2 0];	
%cfg.baselinetype = 'absolute';
cfg.maskstyle    = 'saturation';	
cfg.zlim         = 'maxmin';
cfg.channel      = 'E86';
figure 
ft_singleplotTFR(cfg, FTEEG);
figure 
ft_singleplotTFR(cfg, FTEEG2);

cfg = [];
%cfg.baseline     = [-0.2 0];		
%cfg.baselinetype = 'absolute';
cfg.xlim         = [0 0.4]; 
cfg.ylim         = [15 20];  
cfg.zlim         = 'maxmin';
cfg.marker       = 'on';
figure 
ft_topoplotTFR(cfg, FTEEG2);
figure 
ft_topoplotTFR(cfg, FTEEG);




%% plot using EEGlab functions

%plot FT data
figure
tfdata = permute(inddata,[3 2 1]);
tftopo(tfdata,FTEEG.time*1000,freqs,'chanlocs',EEG.chanlocs,'timefreqs',[-150 0 14 20],'showchan',[],'limits',[-200 400 4 40 -10 10]);

% plot EEGlab data
figure
tfdata = permute(inddata,[3 2 1]);
tftopo(tfdata,times,freqs,'chanlocs',EEG.chanlocs,'timefreqs',[0 400 14 20],'showchan',[],'limits',[-50 400 4 40 -10 10]);

