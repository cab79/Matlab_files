function EEG = acstp_on_epochs(EEG,elec,eventtypes,maskArt,maskTarg,subtractArt,num_trials)

if ~isempty(num_trials); EEG = pop_select(EEG,'trial',1:num_trials);end;

figure
plot(EEG.xmin:0.001:EEG.xmax, mean(EEG.data(elec,:,:),3));title('original');

if ~isempty(maskArt)
    [EEGart, ACSTPstruct, com] = pop_acstp(EEG,'targetevents',eventtypes, 'nontargetevents',{},'timelimits', [EEG.xmin EEG.xmax], 'masktime', maskArt,'subspacedim', [floor(size(EEG.data,1)/2):-1:floor(size(EEG.data,1)/4)],'plotstats','on');
    figure
    plot(EEGart.xmin:0.001:EEGart.xmax, mean(EEGart.data(elec,:,:),3));title('artefact');
end

if subtractArt==1
    EEG.data(:,1:size(EEGart.data,2),:)=EEG.data(:,1:size(EEGart.data,2),:)-EEGart.data;
end

[EEGp, ACSTPstruct, com] = pop_acstp(EEG,'targetevents',eventtypes, 'nontargetevents',{},'timelimits', [EEG.xmin EEG.xmax], 'masktime', maskTarg,'maskelec',elec,'subspacedim', [floor(size(EEG.data,1)/2):-1:floor(size(EEG.data,1)/4)],'plotstats','on');
figure
plot(EEGp.xmin:0.001:EEGp.xmax, mean(EEGp.data(elec,:,:),3));title('corrected');

EEG=EEGp;