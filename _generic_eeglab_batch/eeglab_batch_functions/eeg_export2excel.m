function S=eeg_export2excel(S)

if 0
% SAVE (this part needs to be outside of the loop to save all as one file)
%av_area = permute(av_area,[1 3 2]); % Produces output which is region 1, freq 1; region 1, freq2 etc. Comment out if you want freq1, region 1; freq1, region 2 etc.
av_area = reshape(av_area,[],size(av_area,2)*size(av_area,3)); % columns: freq, regions
av_cell = horzcat(entrainname,num2cell(av_area));
save('averages.mat','av_area','entrainname','av_cell','grand_av_e','grand_av_b','grand_av'); % save filenames and numbers all in one cell array
end