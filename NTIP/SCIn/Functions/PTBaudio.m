function h = PTBaudio(h,opt)


switch opt
    case 'setup'
        
    %configure soundcard
    InitializePsychSound(1);
    h.edevices = PsychPortAudio('GetDevices');
    h.DeviceIndex = [h.edevices.DeviceIndex];
    h.DeviceN = h.DeviceIndex(strcmp({h.edevices.HostAudioAPIName},'ASIO'));
    if isempty(h.DeviceN)
        %outdevices = find(cell2mat({h.edevices.NrOutputChannels})>0 & ~cellfun(@isempty,(strfind({h.edevices.DeviceName},'Speakers'))));
        %h.DeviceN = h.DeviceIndex(outdevices(2));
        h.DeviceN = []; % default sound device
    end
        
    
    
    try
        % (1) =  sound device
        % (2) 1 = sound playback only
        % (3) 1 = default level of latency
        % (4) Requested frequency in samples per second
        % (5) 2 = number of channels
        h.pahandle = PsychPortAudio('Open', h.DeviceN, 1, 1, h.Settings.fs, h.Settings.nrchannels);
    catch
        % Failed. Retry with default frequency as suggested by device:
        fprintf('\nCould not open device at wanted playback frequency of %i Hz. Will retry with device default frequency.\n', freq);
        fprintf('Sound may sound a bit out of tune, ...\n\n');%%
        psychlasterror('reset');
        h.pahandle = PsychPortAudio('Open', h.DeviceN, [], 0, [], h.Settings.nrchannels);
        %pahandle = 1;
    end
end