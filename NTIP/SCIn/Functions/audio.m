function h = audio(h)

if strcmp(h.Settings.stimcontrol,'audioplayer')
    opt = 'audioplayer';
elseif strcmp(h.Settings.stimcontrol,'PsychPortAudio')
    opt = 'PsychPortAudio';
else
    opt = 'audioplayer';
end

% Load audio
t1=GetSecs;
% create index of audio files played, to enable multiple files in
% h.Settings.audiofile
if ~isfield(h,'audiofileindex')
    h.audiofileindex =0;
end
h.audiofileindex = h.audiofileindex + 1;
filename = fullfile(h.d.root,'audio_files',h.Settings.audiofile{h.audiofileindex});
[y,Fs] = audioread(filename);
if isfield(h.Settings,'audiochan')
    ytemp = zeros(size(y,1),h.Settings.nrchannels);
    ytemp(:,h.Settings.audiochan) = y;
    y=ytemp;
end

switch opt
    
    case 'audioplayer'
        h.aud = audioplayer(y, Fs);
        t2=GetSecs;
        disp(['Delay to load audio: ' num2str(t2-t1) ' s'])
        playblocking(h.aud);
    
    case 'PsychPortAudio'
        
        if isfield(h,'pahandle')
            if Fs~=h.Settings.fs || size(y,2)~=h.Settings.nrchannels
                reset=1;
            end
        end
        h = PTBaudio(h,Fs,size(y,2));
        PsychPortAudio('FillBuffer', h.pahandle, y');
        
        startCue = 0;

        % Start audio playback
        waitForDeviceStart = 1;
        PsychPortAudio('Start', h.pahandle, 1, startCue, waitForDeviceStart);
        PsychPortAudio('Stop', h.pahandle, 1);
        
        if reset
            h = PTBaudio(h);
        else
            PsychPortAudio('Close', h.pahandle);
        end
            
end
