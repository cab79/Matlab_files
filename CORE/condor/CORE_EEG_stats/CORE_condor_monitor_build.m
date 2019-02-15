function CORE_condor_monitor_build

nowtime = now;

fin=0;
while fin==0
    pause(5)
    
    if exist('CORE_condor_EEG_job.exe','file')
        s = dir('CORE_condor_EEG_job.exe');
        disp(['file size, bytes: ' num2str(s.bytes)])
        if s.bytes>0 && s.datenum>nowtime
            fin=1;
        end
    end
end
quit