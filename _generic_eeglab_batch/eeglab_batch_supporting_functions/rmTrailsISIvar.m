function EEG = rmTrailsISIvar(EEG,errormarg,timebin,isi,marker)

%dins = [];
%stims = [];
%dintimes = [];
%stimtimes = [];
%for ev = 1:length(EEG.event)
%    if strcmp(EEG.event(ev).type,'DIN2')
%        dins = [dins ev];
%        dintimes = [dintimes EEG.event(ev).init_time];
%    elseif strcmp(EEG.event(ev).type,'STIM')
%        stims = [stims ev];
%        stimtimes = [stimtimes EEG.event(ev).init_time];
%    end
%end
%dinisis = dintimes(2:end)-dintimes(1:end-1);
%stimisis = stimtimes(2:end)-stimtimes(1:end-1);
%dinisis(dinisis<0)=[];
%stimisis(stimisis<0)=[];

% calculate ISIs
etime = nan(1,length(EEG.epoch));
for ep = 1:length(EEG.epoch)
    dinidx = find(strcmp(marker,EEG.epoch(ep).eventtype));
    etime(1,ep) = EEG.epoch(ep).eventinit_time{1,dinidx};
end
dinisis = etime(2:end)-etime(1:end-1);
dinisis_plot = dinisis;
dinisis_plot(dinisis<0)=[];

% plot range of ISI error margins vs trials
pri = 0.001;
tlength = timebin(2)-timebin(1);
pri_range = pri:pri:0.02;
np = [];
for p = pri_range
    np = [np length(find(dinisis_plot<isi-p | dinisis_plot>isi+p))];
end
figure
plot(pri_range,np);

%qd_handle = questdlg('Error margin acceptable at 0.003?', ...
%                 'question', ...
%                 'Yes', 'No', 'Yes');
%switch qd_handle,
%    case 'No',
%      disp('PAUSED - restart if required');
%      pause
%end

%EEGm = EEG;


% create separate EEGs for different error margins - TEST ONLY
%epochs003 = find(dinisis<0.4-0.003 | dinisis>0.4+0.003);
%epochs001 = find(dinisis<0.4-0.001 | dinisis>0.4+0.001);
%[~,ia,ib] = intersect(epochs003,epochs001);
%epochs001003 = epochs001;
%epochs001003(ib) = [];

%EEG3 = pop_select(EEGm,'trial',find(epochs003+1));
%EEG1 = pop_select(EEGm,'trial',find(epochs001003+1));
%selectepochs = ones(1,length(EEG.epoch));
%selectepochs(1) = 0; % remove first trial
%selectepochs(epochs001+1)=0;
%EEG0 = pop_select(EEGm,'trial',find(selectepochs));

%EEG=EEG0;

% remove epochs greater than 'errormarg' of error
rmepochs = find(dinisis<isi-errormarg | dinisis>isi+errormarg);
selectepochs = ones(1,length(EEG.epoch));
selectepochs(1) = 0; % remove first trial
selectepochs(rmepochs+1)=0;
EEG = pop_select(EEG,'trial',find(selectepochs));