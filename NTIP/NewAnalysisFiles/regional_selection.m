%Region/area selection from chanlocs2 file. Doesn't loop, manually had to
%move each file to dname to get them to save in the group. Put all files
%for one participant into file and will save at the end as one file. COM
%values (SeqA-D) are first then UNC values (SeqA-D). 
load chanlocs2

dname = 'Y:\Marie Shorrock\NTIP\Auditory Entrainment Study\Frequency\Region';
entrainfiles = dir(fullfile(dname, '*entrain*'));

cd(dname);

av_area = [];

for f= 1:length(entrainfiles)
    
    entrainname = entrainfiles(f).name;
    basename = strrep(entrainname,'entrain','base');
 
    entrain = load(entrainname); %entrain.f.data
    base = load(basename); % base.fdata
    
    av_e = squeeze(mean(entrain.fdata.powspctrm,1)); %avg powspctrm entrain before correction
    av_b = squeeze(mean(base.fdata.powspctrm,1)); %avg powspctrm base 
    av = av_e./av_b; %divides entrainment powerspectrum by baseline as a correction


    area = [chanlocs2.area];

    uArea = unique(area); %sorts the areas and gets rid of repeats.
    nArea = length(uArea); %find the number of areas

    for r = 1:nArea
           av_area(f,r,:) = mean(av(area,uArea(r),:)); 
    end

    av_area = permute(av_area,[1 3 2]); %this swaps the second and third one
    av_area = reshape(av_area,[],size(av_area,2)*size(av_area,3)); % columns: freq, regions

    
    save NTIP_0001_average
end
