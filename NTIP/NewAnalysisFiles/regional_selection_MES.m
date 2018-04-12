%Region/area selection based on chanlocs2 file. Saves everything at the end
%under the final participant, UNC, SeqD.
dname = 'Y:\Marie Shorrock\NTIP\Auditory Entrainment Study\Frequency';
entrainfiles = dir(fullfile(dname, '*entrain*'));

cd(dname);

av_area = [];
grand_av=[];

for f= 1:length(entrainfiles)
    
    entrainname{f,1} = entrainfiles(f).name; % CAB: made entrainname a cell array of file names, which can be saved at the end.
    basename = strrep(entrainname{f},'entrain','base');
 
    entrain = load(entrainname{f}); %entrain.f.data
    base = load(basename); % base.fdata
    
    av_e = squeeze(mean(entrain.fdata.powspctrm,1)); %avg powspctrm entrain before correction
    av_b = squeeze(mean(base.fdata.powspctrm,1)); %avg powspctrm base 
    av = av_e./av_b; %divides entrainment powerspectrum by baseline as a correction
    grand_av_e = cat(3,grand_av_e,av_e); % CAB: added grand average
    grand_av_b = cat(3,grand_av_b,av_b); % CAB: added grand average
    grand_av = cat(3,grand_av,av); % CAB: added grand average

    area = [chanlocs2.area];

    uArea = unique(area); %sorts the areas and gets rid of repeats.
    nArea = length(uArea); %find the number of areas

    for r = 1:nArea
           av_area(f,r,1:size(av,2)) = squeeze(mean(av(area==uArea(r),:),1))'; % CAB: updated
    end

end

% CAB: This part needs to be outside of the loop to save all as one file:
av_area = permute(av_area,[1 3 2]); %this swaps the second and third one-change if wrong
av_area = reshape(av_area,[],size(av_area,2)*size(av_area,3)); % columns: freq, regions
av_cell = horzcat(entrainname,num2cell(av_area));
grand_av_e = mean(grand_av_e,3); % CAB: added grand average
grand_av_b = mean(grand_av_b,3); % CAB: added grand average
grand_av = mean(grand_av,3); % CAB: added grand average
save('averages.mat','av_area','entrainname','av_cell','grand_av_e','grand_av_b','grand_av'); % CAB: added av_cell - filenames and numbers all in one cell array
