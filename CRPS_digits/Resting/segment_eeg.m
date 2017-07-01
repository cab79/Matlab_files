function EEG = segment_eeg(EEG,dur) % dur = the required duration of the new events in s
name = 'M'; % name of the new events that appears in EEGLAB scroll plot
newEvent = {name 0 dur}; % initiate the new event list (see help pop_importevent for details)
% the following loop will populate the list with required information
for k = 2:round(EEG.xmax/dur) % the loop runs until the end of the loaded dataset
    newEvent = cat(1, newEvent, {name newEvent{k-1,2}+dur dur}); %specs of one event are added on each run of the loop
end
EEG = pop_importevent( EEG, 'event',newEvent,'fields',{'type' 'latency' 'duration'},'timeunit',1,'optimalign','off'); % import the created event list to the loaded dataset