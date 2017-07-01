clear all
loadpaths
files = dir(fullfile(filepath,'*100Hz.Exp3_*subcomp.set'));

for f = 1%:length(files)
    filename = files(f).name;
    [pth basename ext] = fileparts(filename);
    EEG = pop_loadset('filename',filename,'filepath',filepath);
    %[spectra,freqlist,speccomp,contrib,specstd] = pop_spectopo(EEG,1,[],'EEG','freqrange',[0.1 40]);
    [spectra,freqlist,speccomp,contrib,specstd] = pop_spectopo(EEG,1,[],'EEG');
    save(fullfile(filepath, 'spectra', [basename '_spectra.mat']),'spectra','freqlist','speccomp','contrib','specstd');
end


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



