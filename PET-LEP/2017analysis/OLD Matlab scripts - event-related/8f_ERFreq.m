clear all

subjects = {'P1_';'P2_';'P5_';'P6_';'P7_';'P8_';'P14_';'P16_';'P19_';'P20_';'P22_';'P23_';'P24_';'P25_';'P27_';'P30_';'P31_';'P32_';'P33_';'P35_';'S1_';'S2_';'S3_';'S5_';'S6_';'S8_';'S9_';'S10_';'S11_';'S18_';'S20_'};
Nsub = length(subjects);
Ncond = 2;
Nfreq = 3;
ele = {
    [17 18 61]; 
    [13 53 55 47 25 27];
    };
srate = 500; %Hz
timebins = {
    [-3500:(1000/srate):-3000];
    [-2500:(1000/srate):-2000];
    [-500:(1000/srate):0];
    [200:(1000/srate):800];
    }; % in ms relative to pain stimulus
Ntb = length(timebins);
afreqs = cell(Nsub,Nfreq);
results = zeros(Nsub,Ncond,Ntb,Nfreq);

for sub = 1:length(subjects)
    subject = subjects(sub);
    subject = subject{:};
    
    load([char(subjects(sub)) 'total_data_ICA2.mat']);
   
     %%%%%%%%%%%%%%%%%%%%%%%%%%

        % Optional inputs:
        %   'freq'     = [float vector (Hz)] vector of frequencies for topoplot() scalp maps
        %                of power at all channels, or single frequency to plot component 
        %                contributions at a single channel (see also 'plotchan').
        %   'chanlocs' = electrode locations file (format: >> topoplot example)
        %   'limits'   = axis limits [xmin xmax ymin ymax cmin cmax]
        %                To use data limits, omit final values or use nan's
        %                i.e. [-100 900 nan nan -10 10], [-100 900]
        %                Note that default color limits are symmetric around 0 and are
        %                different for each head {defaults: all nans}
        %   'title'    = [quoted string] plot title {default: none}
        %   'freqfac'  = [integer] ntimes to oversample -> frequency resolution {default: 2}
        %   'nfft'     = [integer] value to zero pad data to. Overwrites 'freqfac' above.
        %   'winsize'  = [integer] window size in data points {default: from data}
        %   'overlap'  = [integer] window overlap in data points {default: 0}
        %   'percent'  = [float 0 to 100] percent of the data to sample in computing the 
        %                spectra. Can be used to speed up the computation. {default: 100}.
        %   'freqrange' = [min max] frequency range to plot. Changes x-axis limits. Default is
        %                1 Hz for the min and niquist (srate/2) for the max (if some scalp
        %                maps are plotted, the scalp map at the highest frequency specifies
        %                the maximum).
        %   'reref'    = ['averef'|'off'] convert input data to average reference 
        %                {default: 'off'}
        %   'mapnorm'  = [float vector] if 'data' contain the activity of an independant 
        %                component, this parameter should contain its scalp map. In this case
        %                the spectrum amplitude will be scaled by component RMS scalp power.
        %                Useful for comparing component strengths.
        %   'boundaries' = data point indices of discontinuities in the signal
        %   'plot'     = ['on'|'off'] 'off' -> disable plotting. {default: 'on'}
        %   'rmdc'     =  ['on'|'off'] 'on' -> remove DC. {default: 'off'}  


        %[spectra,freqs,speccomp,contrib,specstd] = spectopo(data, frames, srate, 'key1', 'val1', 'key2', 'val2' ...);
        
        frames = 2750; %frames per epoch {0 -> data length}
        freqrange = [6 14];
        midfreq = 10;
        freqstep = 0.5;
        
        [spectra,freqs,speccomp,contrib,specstd] = spectopo(total_data_ICA2(:,:,:), frames, srate,'winsize',1000,'overlap',500,'freqrange',freqrange);
        close all;
        ele_avg = mean(spectra,1);
            
        %i_max_a = 1;
        fstep = 0;
        max_a=[];
        while isempty(max_a)
            fstep=fstep+1;
            a = find(freqs>(midfreq-freqstep*fstep) & freqs<(midfreq+freqstep*fstep));
            sub_a = find(freqs<=(midfreq-freqstep*fstep));
            [max_a i_max_a] = findpeaks(ele_avg(:,a),'SORTSTR','descend');
     
            if ~isempty(max_a)
                if max(ele_avg(:,a)) ~= max_a(1) 
                    max_a=[];
                end
            end
        end;
        
        i_max_a = i_max_a(1);
        max_index = length(sub_a)+i_max_a;
        IAF = freqs(max_index);
        
        if IAF<freqrange(1) || IAF>freqrange(2)
            IAF = midfreq;
        end
        
        plot(freqs(a),ele_avg(:,a)); hold on; scatter(IAF,max_a(1),'r');title(char(subjects(sub)));
        pause(1);
        a1 = find(freqs>(IAF-4) & freqs<(IAF-2));
        a2 = find(freqs>(IAF-2) & freqs<(IAF));
        a3 = find(freqs>(IAF) & freqs<(IAF+2));
        afreqs{sub,1} = a1;
        afreqs{sub,2} = a2;
        afreqs{sub,3} = a3;
end

for sub = 1:length(subjects)
    subject = subjects(sub);
    subject = subject{:};
    
    load([char(subjects(sub)) 'total_data_ICA2.mat']);
    load([subject 'data_matrix_dim2.mat']);
    eval(['data = total_data_ICA2;']);

    if x(3)+x(4) ~= size(data,3); errordlg('trials per condition do not match data size');end
    t1 = 1:x(3);
    t2 = x(3)+1:x(3)+x(4);
    data1 = data(:,:,t1);
    data2 = data(:,:,t2);

    clear data t1 t2 x
    
    for dt = 1:2

        eval(['data = data' num2str(dt) ';' ]);

        for tb = 1:length(timebins)
            dtime = ((timebins{tb}+4000)/(1000/srate));
            [spectra,freqs,speccomp,contrib,specstd] = spectopo(data(:,dtime,:), length(dtime), srate,'freqrange',freqrange,'plot','off');    
            for e = 1:length(ele)
                avg_spect = mean(spectra(ele{e},:),1);
                results(sub,dt,tb,1,e) = 10.^(mean(avg_spect(:,afreqs{1}),2)/10);
                results(sub,dt,tb,2,e) = 10.^(mean(avg_spect(:,afreqs{2}),2)/10);
                results(sub,dt,tb,3,e) = 10.^(mean(avg_spect(:,afreqs{3}),2)/10);
            end
            %pause;
        end

        clear data
    end

    clear a1 a2 a3;
end

results_2d = reshape(results,length(subjects),2*length(timebins)*3*length(ele));
save alpha_power_f3_tb4_ele2.mat results_2d

