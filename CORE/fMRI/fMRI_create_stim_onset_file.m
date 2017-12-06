%% settings
studyname = 'CORE';
subnum = 'pilot';
dname = 'C:\Matlab_files\CORE\fMRI\SCIn_outputs';
seq_fname = 'Sequence_pilot_Sequence_CORE_fMRI_cont_8block_OptionfMRI_run_startblock1_20171024T104134.mat';
out_fname = 'Output_pilot_Sequence_CORE_fMRI_cont_8block_OptionfMRI_run_startblock1_20171024T104134.mat';

% 'expected' is based on stimulus design; 
% 'recorded' is that recorded by the software - older data has a bug meaning this is not accurate
onset_marker = 'expected';
correct_missing = 1; % if there are missing stim onsets, infer what they should be. Number here is estimated ISI. Only works if ISI is fixed throughout expt.

% event condition names, "condition" numbers, and "signal" numbers of interest
events = {
    'tact_slow_odd',1,[1 2]
    'aud_slow_odd',3,[3 4]
    'tact_fast_odd',5,[1 2]
    'aud_fast_odd',7,[3 4]
    'tact_slow_odd_chan1',1,[1]
    'aud_slow_odd_chan1',3,[3]
    'tact_fast_odd_chan1',5,[1]
    'aud_fast_odd_chan1',7,[3]
    'tact_slow_odd_chan2',1,[2]
    'aud_slow_odd_chan2',3,[4]
    'tact_fast_odd_chan2',5,[2]
    'aud_fast_odd_chan2',7,[4]
    };

% block names and associated condition numbers
blocks = {
    'tact_slow',[1 2],[1 2]
    'aud_slow',[3 4],[3 4]
    'tact_fast',[5 6],[1 2]
    'aud_fast',[7 8],[3 4]
    %'tact_slow_chan1',[1 2],[1]
    %'aud_slow_chan1',[3 4],[3]
    %'tact_fast_chan1',[5 6],[1]
    %'aud_fast_chan1',[7 8],[3]
    %'tact_slow_chan2',[1 2],[2]
    %'aud_slow_chan2',[3 4],[4]
    %'tact_fast_chan2',[5 6],[2]
    %'aud_fast_chan2',[7 8],[4]
    };
blockmarker = 0;

%% RUN
dname = fullfile(dname,[studyname subnum]);
cd(dname)
load(seq_fname)
load(out_fname)

if length(out.trigtime)<4
    error('should be 4 scanner triggers recorded')
end

scanstart = out.trigtime(1);
rectime = out.stimtime;
nostimrec = cellfun(@isempty,rectime);
stimrec = ~nostimrec;
rectime(nostimrec)={NaN};
rectime = cell2mat(rectime);

% fix expisi if too small
expisi = out.expisi;
expisi(1,end:length(stimrec)-1) = {NaN};
expisi = cell2mat(expisi);

switch onset_marker
    case 'recorded'
        stimtime = rectime - scanstart;
        
    case 'expected'
        if correct_missing
            expisi(isnan(expisi)) = correct_missing;
        end
        stimtime = [0 cumsum(expisi)] + rectime(1)-scanstart;
end

for i = 1:size(events,1) 
    conds = ismember(seq.condnum,events{i,2}) & ismember(seq.signal,events{i,3});
    onsets = stimtime(conds)';
    
    %save([events{i,1} '.txt'],'onsets','-ASCII')
    fid = fopen(fullfile(dname,[studyname subnum '_' events{i,1} '_event.txt']), 'wt');
    fprintf(fid, '%.2f\n', onsets);
    fclose(fid);
end

for i = 1:size(blocks,1)
    conds = ismember(seq.condnum,blocks{i,2}) & ismember(seq.signal,blocks{i,3});
    st=find(diff(conds)==1); % don't add 1 because we will start from the first '0' condition
    en=find(diff(conds)==-1)+1;
    
    % add on last trial to en if needed
    if length(en)==length(st)-1
        en = [en length(conds)];
    end
    
    onsets = stimtime(st)';
    durations = stimtime(en)'-stimtime(st)';
    
    fid = fopen(fullfile(dname,[studyname subnum '_' blocks{i,1} '_block.txt']), 'wt');
    fprintf(fid, '%.2f\n', onsets);
    fclose(fid);
    
    fid = fopen(fullfile(dname,[studyname subnum '_' blocks{i,1} '_block_dur.txt']), 'wt');
    fprintf(fid, '%.2f\n', durations);
    fclose(fid);
end