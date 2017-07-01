% from a design template, reverses the block order and saves a new version

clear all
dname=('C:\Users\cab79\Google Drive\2. Cambridge University\programs\Digitimer programs via labjack\Design templates');
fname = 'TSOT_design_template_part4_twohands_balanced';

load(fullfile(dname,fname));

[hand,dc,cp,bi,blockii,btypes] = blocktype(dname,fname);

bii = fliplr([blockii size(design,2)+1]);
rev=[];
for b = 1:length(bi)
    rev = [rev design(1:3,bii(b+1):bii(b)-1)];
end
rev(4,:) = design(4,:);

design=rev;
s_name=fullfile(dname,[fname '_reversed']);
save(s_name,'design');