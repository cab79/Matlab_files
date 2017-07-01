%Aim: based on a list of behavioural events (markers 1 - 6), name the
%trials according to the index of those behavioural events. 
%Produces: columns where row index indicates marker number in EEG, and the
%value in that row indicates the index from the list of condition numbers
%in the XX_conds.mat file. That list of condition numbers is the trials
%that are fed into the HGF model.

clear all
filepath = 'C:\Data\Expectancy Study\Preprocessed';
condpath = 'M:\Matlab\ExpStudy\Behaviour\Results\cond numbers';
rfile = 'ExpStudy_behav_idx.xlsx';
cd(filepath);
files = dir('*orig.set');
T = table();

for f = 1:length(files)
    orig_file = files(f).name;
    [pth nme ext] = fileparts(orig_file); 
    C = strsplit(nme,'_');
    if strcmp(C{1},'F4'); continue; end;
    EEG = pop_loadset('filename',orig_file,'filepath',filepath);
    load(fullfile(condpath,[C{1} '_conds.mat']));
    markers = cellfun(@str2double,{EEG.event.type})';
    
    del = 0;
    for j = 1:size(EEG.event,2)
        if ~(strcmp(EEG.event(j-del).type, '1') || strcmp(EEG.event(j-del).type, '2') || strcmp(EEG.event(j-del).type, '3') || strcmp(EEG.event(j-del).type, '4') || strcmp(EEG.event(j-del).type, '5') || strcmp(EEG.event(j-del).type, '6'))
            EEG.event(j-del) = [];
            del = del + 1;
        end
    end
    
    % for each marker, identify the corresponding index of cond
    cond_idx = NaN(length(markers),2); 
    for m = 1:length(markers) % for each marker...
        if isnan(cond_idx(m)) % ...we need to make sure there is an associated value in cond_idx
            for n = 1:length(markers)-m % if not, try every possible length of a sequence (mseq), starting from the smallest...
                mseq = markers(m:m+n);
                match=0;
                for c =1:length(cond)-length(mseq)+1 % ...and considering all possible values of cond...
                    if sum(abs(mseq-cond(c:c+length(mseq)-1)))==0 % ...find the number of times this sequence occurs.
                        match = match+1;
                        cind = c;
                    end
                end
                if match>1 && length(mseq)<10
                    continue
                elseif match==0
                    break
                else
                    %cond_idx(m,1) = cind;
                    cond_idx(m:m+length(mseq)-1,1) = cind:cind+length(mseq)-1;
                    cond_idx(m,2) = n+1;
                end
            end
        end
    end
    
    % remove short sequences placed incorrectly
    isflags=1;
    flagcount=0;
    while isflags==1
        flagcount=flagcount+1
        flag_idx=zeros(size(cond_idx,1),1);
        condii = [1:size(cond_idx,1)]';
        nonan_cond_idx = cond_idx(~isnan(cond_idx(:,1)),:);
        nonan_condii = condii(~isnan(cond_idx(:,1)),:);
        for m = 2:size(nonan_cond_idx,1)-1 % identify where the sequence decreases when it should be increasing
            if nonan_cond_idx(m-1,1)>nonan_cond_idx(m,1)
                flag_idx(nonan_condii(m),1)=1;
            end
            %if ~isnan(cond_idx(m+1,1)) && ~isnan(cond_idx(m,1))
            %    if cond_idx(m+1,1)<cond_idx(m,1)
            %        flag_idx(m,1)=1;
            %    end
            %end
        end
        flag_idx=find(flag_idx);
        if isempty(flag_idx); isflags=0;end;
        seqst = cond_idx(~isnan(cond_idx(:,2)),1); % find the start of all sequences
        seqst_len = cond_idx(all(~isnan(cond_idx(:,:)),2),2); 
        seqst_idx = find(all(~isnan(cond_idx(:,:)),2)); % find the start of all sequences
        for fg = 1:length(flag_idx)
            fi = find(seqst_idx==flag_idx(fg));
            si = fi-1:fi+1;
            si = intersect(si,1:length(seqst_len));
            [si_min,si_min_idx] = min(seqst_len(si));
            repst = seqst_idx(si_min_idx+min(si)-1);
            if si_min<10
                cond_idx(repst:repst+si_min-1) = NaN;
            else
                flag_idx(fg)=[];
            end
        end
        if isempty(flag_idx); isflags=0;end;
    end
    
    tab = [cond_idx(1:min(720,size(cond_idx,1)),1); NaN(720-size(cond_idx,1),1)];
    
    T.(C{1})=tab;
    
end

writetable(T,rfile);