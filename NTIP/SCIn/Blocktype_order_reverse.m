% from a design template, reverses the block order and saves a new version

clear all
dname=('C:\Matlab_files\NTIP\SCIn\Sequences');
fname = 'Sequence_CORE_fMRI_cont_8block_OptionfMRI_run';

load(fullfile(dname,fname));

[mode,cp,bi,blockii,btypes] = blocktype_CORE_fMRI(dname,fname);

bii = fliplr([blockii size(seq.condnum,2)+1]);
rev.signal=[];
rev.condnum=[];
rev.changedist=[];
rev.blocks=[];
rev.RP=[];
for b = 1:length(bi)
    rev.signal = [rev.signal seq.signal(bii(b+1):bii(b)-1)];
    rev.condnum = [rev.condnum seq.condnum(bii(b+1):bii(b)-1)];
    rev.changedist = [rev.changedist seq.changedist(bii(b+1):bii(b)-1)];
    rev.blocks = [rev.blocks seq.blocks(bii(b+1):bii(b)-1)];
    rev.RP = [rev.RP seq.RP(bii(b+1):bii(b)-1)];
end

seq=rev;
s_name=fullfile(dname,[fname '_reversed']);
save(s_name,'seq','settings');