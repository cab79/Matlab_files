clear all

subjects = {'P1_';'P2_';'P5_';'P6_';'P7_';'P8_';'P15_';'P16_';'S1_';'S2_';'S3_';'S5_';'S6_';'S8_';'S9_';'S10_';'S11_';'S18_';'S20_';'S4_';'P14_';'P20_';'P24_';'P30_';'P32_';'P33_';'P19_';'P22_';'P23_';'P27_';'P31_'};
ele = 18;;
srate = 500; %Hz
tlimits = [-1000 1500]; % in ms relative to pain stimulus
baseline_end = [0];
cycles = [3 0.5];
total_data = [];

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
        total_data = cat(3,total_data,data);

    end
end

        frames=((tlimits(2)-tlimits(1))/2)+1;
        ntrials=size(total_data,3);
        inputdata = reshape(mean(total_data(ele,(tlimits(1)+4000)/2:(tlimits(2)+4000)/2,:),1),1,frames*ntrials);
        %%%%%%%%%%%%%%%%%%%%%%%%%%
        %Required inputs:     
        %       data        = Single-channel data vector (1,frames*ntrials) (required)
        %       frames      = Frames per trial                     {def|[]: datalength}
        %       tlimits     = [mintime maxtime] (ms) Epoch time limits 
        %                      {def|[]: from frames,srate}
        %       srate       = data sampling rate (Hz)                  {def:250}
        %       cycles      = If 0 -> Use FFTs (with constant window length) {0 = FFT}
        %                     If >0 -> Number of cycles in each analysis wavelet 
        %                     If [wavecycles factor] -> wavelet cycles increase with 
        %                     frequency  beginning at wavecyles (0<factor<1; factor=1 
        %                     -> no increase, standard wavelets; factor=0 -> fixed epoch 
        %                     length, as in FFT.  Else, 'mtaper' -> multitaper decomp. 
        %
        %    Optional Inter-Irial Coherence (ITC) type:
        %       'type'      = ['coher'|'phasecoher'] Compute either linear coherence 
        %                      ('coher') or phase coherence ('phasecoher') also known
        %                      as the phase coupling factor           {'phasecoher'}.
        %    Optional detrending:
        %       'detret'    = ['on'|'off'], Detrend data in time.               {'off'}
        %       'detrep'    = ['on'|'off'], Detrend data across trials          {'off'}
        %
        %    Optional FFT/DFT parameters:
        %       'winsize'   = If cycles==0: data subwindow length (fastest, 2^n<frames);
        %                     If cycles >0: *longest* window length to use. This
        %                      determines the lowest output frequency       {~frames/8}
        %       'timesout'  = Number of output times (int<frames-winframes)       {200}
        %       'padratio'  = FFT-length/winframes (2^k)                            {2}
        %                      Multiplies the number of output frequencies by
        %                      dividing their spacing. When cycles==0, frequency
        %                      spacing is (low_freq/padratio).
        %       'maxfreq'   = Maximum frequency (Hz) to plot (& to output if cycles>0) 
        %                      If cycles==0, all FFT frequencies are output.      {50}
        %       'baseline'  = Spectral baseline window center end-time (in ms).    {0}
        %       'powbase'   = Baseline spectrum (power, not dB) to normalize the data. 
        %                      {def|NaN->from data}
        %
        %    Optional multitaper parameters:
        %       'mtaper'    = If [N W], performs multitaper decomposition. 
        %                      (N is the time resolution and W the frequency resolution; 
        %                      maximum taper number is 2NW-1). Overwrites 'winsize' and 
        %                      'padratio'. 
        %                     If [N W K], uses K Slepian tapers (if possible).
        %                      Phase is calculated using standard methods.
        %                      The use of mutitaper with wavelets (cycles>0) is not 
        %                      recommended (as multiwavelets are not implemented). 
        %                      Uses Matlab functions DPSS, PMTM.      {no multitaper}
        %
        %    Optional bootstrap parameters:
        %       'alpha'     = If non-0, compute two-tailed bootstrap significance prob. 
        %                     level. Show non-signif. output values in green       {0}
        %       'naccu'     = Number of bootstrap replications to accumulate       {200}
        %       'baseboot'  = Bootstrap baseline to subtract (1 -> use 'baseline'(above)
        %                                                     0 -> use whole trial) {1}
        %    Optional scalp map:
        %       'topovec'   = Scalp topography (map) to plot                     {none}
        %       'elocs'     = Electrode location file for scalp map   
        %                     File should be ascii in format of  >> topoplot example   
        %                     May also be an EEG.chanlocs struct. 
        %                     {default: file named in icadefs.m}
        %    Optional plotting parameters:
        %       'hzdir'     = ['up'|'down'] Direction of the frequency axes; reads default
        %                     from icadefs.m                                     {'up'}
        %       'plotersp'  = ['on'|'off'] Plot power spectral perturbations     {'on'} 
        %       'plotitc'   = ['on'|'off'] Plot inter trial coherence            {'on'}
        %       'plotphase' = ['on'|'off'] Plot sign of the phase in the ITC panel, i.e.
        %                     green->red, pos.-phase ITC, green->blue, neg.-phase ITC {'on'}
        %       'erspmax'   = [real dB] set the ERSP max. for the scale (min= -max){auto}
        %       'itcmax'    = [real<=1] set the ITC maximum for the scale          {auto}
        %       'title'     = Optional figure title                                {none}
        %       'marktimes' = Non-0 times to mark with a dotted vertical line (ms) {none}
        %       'linewidth' = Line width for 'marktimes' traces (thick=2, thin=1)  {2}
        %       'pboot'     = Bootstrap power limits (e.g., from timef())    {from data}
        %       'rboot'     = Bootstrap ITC limits (e.g., from timef())      {from data}
        %       'axesfont'  = Axes text font size                                  {10}
        %       'titlefont' = Title text font size                                 {8}
        %       'vert'      = [times_vector] -> plot vertical dashed lines at given ms.
        %       'verbose'   = ['on'|'off'] print text                              {'on'}
        %
        %    Outputs: 
        %            ersp   = Matrix (nfreqs,timesout) of log spectral diffs. from 
        %                     baseline (in dB).  NB: Not masked for significance. 
        %                     Must do this using erspboot
        %            itc    = Matrix of inter-trial coherencies (nfreqs,timesout) 
        %                     (range: [0 1]) NB: Not masked for significance. 
        %                     Must do this using itcboot
        %          powbase  = Baseline power spectrum (NOT in dB, used to norm. the ERSP)
        %            times  = Vector of output times (sub-window centers) (in ms)
        %            freqs  = Vector of frequency bin centers (in Hz)
        %         erspboot  = Matrix (2,nfreqs) of [lower;upper] ERSP significance diffs
        %          itcboot  = Matrix (2,nfreqs) of [lower;upper] ITC thresholds (not diffs)
        %          itcphase = Matrix (nfreqs,timesout) of ITC phase (in radians)

        

        [ersp,itc,powbase,times,freqs,erspboot] = timef(inputdata,frames,tlimits,srate,cycles,'baseline',baseline_end,'alpha',0.05,'maxfreq',20);
        %[ersp,itc,powbase,times,freqs,erspboot] = timef(inputdata,frames,tlimits,srate,cycles,'baseline',baseline_end,'maxfreq',100);
