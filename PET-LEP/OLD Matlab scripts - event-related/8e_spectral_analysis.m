clear all

subjects = {'S1_';'S4_';'S5_';'S8_';'S9_';'S10_';'S11_'};
scan_ord = [2 1 2 1 2 2 1]; % which scan was painful
scans = {'pain','nonpain'};
sessions = {'2','3'};

all_spect = zeros(62,length(subjects),length(sessions),2);
delta = zeros(62,length(subjects),length(sessions),2);
theta = zeros(62,length(subjects),length(sessions),2);
alpha = zeros(62,length(subjects),length(sessions),2);
beta = zeros(62,length(subjects),length(sessions),2);
gamma1 = zeros(62,length(subjects),length(sessions),2);
gamma2 = zeros(62,length(subjects),length(sessions),2);

%load alpha
%load delta
%load theta
%load beta
%load gamma1
%load gamma2

for sess = 1:length(sessions)
    

for sub = 1:length(subjects)
    subject = subjects(sub);
    subject = subject{:};
    
    if sess == scan_ord(sub)
        session = sessions(1);
        session = session{:};
    elseif sess ~= scan_ord(sub)
        session = sessions(2);
        session = session{:};
    end

    
%load([char(subjects(sub)) 'total_data_ICA' session '_acc.mat']);
%load([subject 'data_dim_acc' session]);
%eval(['data = total_data_ICA' session ';']);

load([char(subjects(sub)) 'total_data_ICA' session '.mat']);
load([subject 'data_matrix_dim_' session]);
eval(['data = total_data_ICA' session ';']);
x(4) = 120;
x(5) = 120;

if x(4)+x(5) ~= size(data,3); errordlg('trials per condition do not match data size');end
t1 = 1:x(4);
t2 = x(4)+1:x(4)+x(5);
data1 = data(:,:,t1);
data2 = data(:,:,t2);

clear data t1 t2 t3 x

for dt = 1:2
    
eval(['data = data' num2str(dt) ';' ]);

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

frames = 2500; %frames per epoch {0 -> data length}
srate  = 500; %sampling rate per channel (Hz)
[spectra,freqs,speccomp,contrib,specstd] = spectopo(data, frames, srate,'winsize',1000,'overlap',500);

s = find(freqs>0.5 & freqs<90);
all_spect(:,sub,sess,dt) = 10.^mean(spectra(:,s),2);

d = find(freqs>0.5 & freqs<4);
delta(:,sub,sess,dt) = log10((10.^mean(spectra(:,d),2))*100/mean(all_spect(:,sub,sess,dt),1));
t = find(freqs>4 & freqs<7.5);
theta(:,sub,sess,dt) = log10((10.^mean(spectra(:,t),2))*100/mean(all_spect(:,sub,sess,dt),1));
a = find(freqs>8 & freqs<13);
alpha(:,sub,sess,dt) = log10((10.^mean(spectra(:,a),2))*100/mean(all_spect(:,sub,sess,dt),1));
b = find(freqs>14 & freqs<30);
beta(:,sub,sess,dt) = log10((10.^mean(spectra(:,b),2))*100/mean(all_spect(:,sub,sess,dt),1));
g1 = find(freqs>30 & freqs<60);
gamma1(:,sub,sess,dt) = log10((10.^mean(spectra(:,g1),2))*100/mean(all_spect(:,sub,sess,dt),1));
g2 = find(freqs>60 & freqs<90);
gamma2(:,sub,sess,dt) = log10((10.^mean(spectra(:,g2),2))*100/mean(all_spect(:,sub,sess,dt),1));

clear data d t a b g1 g2
end
end
end

eval(['save gamma1.mat gamma1']);
eval(['save gamma2.mat gamma2']);
eval(['save beta.mat beta']);
eval(['save alpha.mat alpha']);
eval(['save theta.mat theta']);
eval(['save delta.mat delta']);
