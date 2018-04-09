%Region/area selection based on chanlocs2 file. Saves everything at the end
%under the final participant, UNC, SeqD.
dname = 'Y:\Marie Shorrock\NTIP\Auditory Entrainment Study\Frequency';
entrainfiles = dir(fullfile(dname, '*entrain*'));

cd(dname);

av_area = [];


for f= 1:length(entrainfiles)
    
    entrainname{f} = entrainfiles(f).name; % CAB: made entrainname a cell array of file names, which can be saved at the end.
    basename = strrep(entrainname{f},'entrain','base');
 
    entrain = load(entrainname{f}); %entrain.f.data
    base = load(basename); % base.fdata
    
    av_e = squeeze(mean(entrain.fdata.powspctrm,1)); %avg powspctrm entrain before correction
    av_b = squeeze(mean(base.fdata.powspctrm,1)); %avg powspctrm base 
    av = av_e./av_b; %divides entrainment powerspectrum by baseline as a correction


    area = [chanlocs2.area];

    uArea = unique(area); %sorts the areas and gets rid of repeats.
    nArea = length(uArea); %find the number of areas

    for r = 1:nArea
           av_area(f,r,:) = mean(av(area==uArea(r),:),1); % CAB: updated
    end

    av_area = permute(av_area,[1 3 2]); %this swaps the second and third one-change if wrong
    av_area = reshape(av_area,[],size(av_area,2)*size(av_area,3)); % columns: freq, regions
    
    [pth nme ext] = fileparts(entrainname{f}); 
    sname_ext = ['average.mat'];
    sname = [nme '_' sname_ext];
    save(fullfile(sname),'av_area','entrainname'); % CAB: added entrainname
end



