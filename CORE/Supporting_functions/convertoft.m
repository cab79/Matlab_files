%% convert from EEGLAB to FIELDTRIP structure
function EEG = convertoft(EEG)

EEG = eeglab2fieldtrip(EEG,'preprocessing','none');
EEG.elec.pnt = [-EEG.elec.pnt(:,2) EEG.elec.pnt(:,1) EEG.elec.pnt(:,3)];
EEG.elec.elecpos = EEG.elec.pnt;
EEG.elec.chanpos = EEG.elec.elecpos;
EEG.elec.unit = 'cm';