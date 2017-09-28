load('C:\Data\Matlab\Matlab_files\NTIP\SCIn\wavetest.mat')

h.Settings.fs=96000;
h.Settings.nrchannels=2;
%h = PTBaudio(h,[],[],1);
h = PTBaudio(h);
buff = PsychPortAudio('CreateBuffer', h.pahandle, wavetest);
%PsychPortAudio('Start', h.pahandle, 0, 0, 1);

%h.pa1 = PsychPortAudio('OpenSlave', h.pahandle, 1);
%h.pa2 = PsychPortAudio('OpenSlave', h.pahandle, 1);

% Create audio buffers prefilled with the 3 sounds:
% eventually need to fill these buffers on trial-1
%pabuffer1 = PsychPortAudio('CreateBuffer', [], wavetest);
%pabuffer2 = PsychPortAudio('CreateBuffer', [], wavetest);

% Demo of sound schedules: Build a sound schedule (a "playlist") that will
% play the sound buffers.
% The schedule shall contain exactly 2 slots:
PsychPortAudio('UseSchedule', h.pahandle, 1, 750);
PsychPortAudio('AddToSchedule', h.pahandle, buff);
st = PsychPortAudio('Start', h.pahandle, 0, 0, waitForDeviceStart);

keyIsDown=0;
delay = [];
i=0;
startCue = 0;
while i<8
    %pause(0.001)
    %keyIsDown = KbCheck;
    i = i+1;
    
    waitForDeviceStart = 1;
    commandtime = GetSecs;
    try
        commanddelay = commandtime-t2;
    end
    try
        old_st=st;
    catch
        old_st=0;
    end
    % Add next trial to sound buffer
    buff = PsychPortAudio('CreateBuffer', h.pahandle, wavetest);
    PsychPortAudio('AddToSchedule', h.pahandle, buff);
    %st = PsychPortAudio('Start', h.pahandle, 0, 0, waitForDeviceStart);
    % stop playback at end of each trial
    %PsychPortAudio('Start', h.pahandle, 1, startCue, waitForDeviceStart);
    %status = PsychPortAudio('GetStatus', h.pahandle);
    
    % status.starttime = time the playback started, but not when hit speakers
    % status.CurrentStreamTime = Estimate of when the most recently submitted sample will hit the speaker.
    
    %st = status.StartTime;
    %cst = status.CurrentStreamTime;
    %ps = status.PositionSecs;
    %pl = status.PredictedLatency; % of hardware
    
    isi = st-old_st
    delay(i,1) = st-commandtime; % delay in computer starting playback; not including latency to hit speakers
    meandelay = mean(delay);
    
    % estimates of required delay
    %est1 = ps;
    %est2 = meandelay;
    %est3 = meandelay+pl;
    
    t2 = GetSecs;
    %playtime = t2-st; % shorter than 0.8 due to start delay; should get even shorter with startCue>0
    %playtimefromcommand = t2-commandtime; % stops it 0.8s after command time
    %estendtime = st+size(wavetest,2)/status.SampleRate;
    %actendtime = st+playtime;
    %startCue = actendtime-est3; % can't subtract from the current time! Need to stop playing sooner...
end
PsychPortAudio('Stop', h.pahandle, 1);
PsychPortAudio('Close', h.pahandle)