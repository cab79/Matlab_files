function h = PTBaudio(h,varargin)
global pah
Ncset=0;Fsset=0;master=0;
% varargin = sampling rate that is different from settings
if ~isempty(varargin)
    if length(varargin)>2
        master = varargin{3};
    end
    if length(varargin)>1
        h.Nc = varargin{2};
        Ncset=1;
    end
    if length(varargin)>0
        h.Fs = varargin{1};
        Fsset=1;
    end
end

if ~isfield(h,'Fs') || isempty(h.Fs) || ~Fsset
    h.Fs = h.Settings.fs;
end
if ~isfield(h,'Nc') || isempty(h.Nc) || ~Ncset
    h.Nc = h.Settings.nrchannels;
end

%configure soundcard
if ~isfield(h,'pahandle')
    setpref('dsp','portaudioHostApi',3) 
    InitializePsychSound(1);
    h.edevices = PsychPortAudio('GetDevices');
    h.DeviceIndex = [h.edevices.DeviceIndex];
    h.DeviceN = h.DeviceIndex(strcmp({h.edevices.HostAudioAPIName},'ASIO'));
    if isempty(h.DeviceN)
        %outdevices = find(cell2mat({h.edevices.NrOutputChannels})>0 & ~cellfun(@isempty,(strfind({h.edevices.DeviceName},'Speakers'))));
        %h.DeviceN = h.DeviceIndex(outdevices(2));
        h.DeviceN = []; % default sound device
    end
else
    try
        PsychPortAudio('Close', h.pahandle);
    end
    try 
        PsychPortAudio('Close',pah);
    end
end

try
    % (1) =  sound device
    % (2) 1 = sound playback only
    % (3) 1 = default level of latency
    % (4) Requested frequency in samples per second
    % (5) 2 = number of channels
    if master
        h.pahandle = PsychPortAudio('Open', h.DeviceN, 1+8, 1, h.Fs, h.Nc);
    else
        h.pahandle = PsychPortAudio('Open', h.DeviceN, 1, 1, h.Fs, h.Nc);
    end
    pah = h.pahandle;
catch
    try
        PsychPortAudio('Close', h.pahandle);
    end
    try 
        PsychPortAudio('Close',pah);
    end
    try 
        PsychPortAudio('Close');
    end
    h = PTBaudio(h);
    % Failed. Retry with default frequency as suggested by device:
    %fprintf('\nCould not open device at wanted playback frequency of %i Hz. Will retry with device default sampling rate.\n', h.Fs);
    %fprintf('Sound may sound a bit out of tune, ...\n\n');%%
    %psychlasterror('reset');
    %h.pahandle = PsychPortAudio('Open', h.DeviceN, [], 0, [], h.Nc);
    %pahandle = 1;
end