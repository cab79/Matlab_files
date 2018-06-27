
%% OLD VERSION!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
clear all
filepath = 'C:\Data\CORE\Preprocessed_100Hz';
cd(filepath);
files = dir('*_1st_ICA.set');
load('C:\Data\Matlab\Matlab_files\CORE\Supporting_functions\chanlocs.mat');

ALLEEG=struct;
ALLEEG_save = 1; % 1= save multiple ERPs from one file; 2 = save one ERP from multiple files

basebins = [-0.2 0; % for epoching, TSOT(2)
            -0.05 0]; % for epoching, TSOT(4)
        
files_ana = [5 6 83:92]% 1:length(files);
for f = files_ana
    [pth nme ext] = fileparts(files(f).name); 
    C = strsplit(nme,'_');
    EEG = pop_loadset('filename',files(f).name,'filepath',filepath);
    EEG = pop_reref( EEG, []);
    EEG = pop_subcomp( EEG, [], 0);
    EEG = eeg_interp(EEG,eeg_mergelocs(chanlocs),'spherical');
    for i = 1:EEG.trials, EEG.data(:,:,i) = detrend(EEG.data(:,:,i)')'; end;
    
    
    
    EEG = FTrejman(EEG,[0 0]);
    EEG = eeg_interp(EEG,eeg_mergelocs(chanlocs),'spherical');
    %EEG = pop_autorej(EEG, 'nogui','on','threshold',1000,'startprob',15);
    sname = [C{1} '_' C{2} '_cleaned.set'];
    EEG = pop_saveset(EEG,'filename',sname,'filepath',filepath); 
    
end

for f = files_ana
    [pth nme ext] = fileparts(files(f).name); 
    C = strsplit(nme,'_');
    EEG = pop_loadset('filename',[C{1} '_' C{2} '_cleaned.set'],'filepath',filepath);
    
    if strcmp(C{2},'2')
        basebin = basebins(1,:);
    elseif strcmp(C{2},'4')
        basebin = basebins(2,:);
    end
    
    EEG = pop_rmbase(EEG, basebin*1000);
    
    [conds, tnums, fnums, bnums] = get_markers(EEG);
    
    %isis = etime(2:end)-etime(1:end-1);
    %isis(isis<0)=[];
    %plot(isis)
    
    % IF digits analysis is needed and there is no digit (FNUM) along with
    % CNUM:
    %bl=[];
    %for i = 2:length(tnums)
    %    if tnums(i)<tnums(i-1)
    %        bl = [bl i];
    %    end
    %end
    %fnums=[];
    %dt.fnum1 = dt.design(1,:);
    %dt.fnum2 = dt.design(1,[find(dt.design(4,:)==2) find(dt.design(4,:)==3)]);
    %dt.fnum3 = dt.design(1,find(dt.design(4,:)==3));
    %for i = 1:length(tnums(1:bl(1)-1))
    %    fnums = [fnums dt.fnum1(tnums(i))];
    %end
    %for i = 1:length(tnums(bl(1):bl(2)-1))
    %    fnums = [fnums dt.fnum2(tnums(bl(1)+i-1))];
    %end
    %for i = 1:length(tnums(bl(2):end))
    %    fnums = [fnums dt.fnum3(tnums(bl(2)+i-1))];
    %end
    
    % Find indices of all "no change" conditions, if these needs selecting
    % later
    no_change_cond_num = [3 4 7 8 11 12 15 16 19 20 23 24];
    %no_change_cond_num = [3 4 11 12 19 20]; %left only
    
    
    % create an index of conds of the "no change" trials
    no_change_conds=[];
    for i = 1:length(no_change_cond_num)
        no_change_conds = [no_change_conds find(conds==no_change_cond_num(i))];
    end
    no_change_conds = sort(no_change_conds);
    
    % select event type
    ep_mark=conds;
    %ep_mark=fnums;
    select_nochange=0;
    
    EEG_ALL = EEG;
    if ALLEEG_save==1
        ALLEEG=struct;
    end
    ctype = sort(unique(ep_mark));
    ctype(ctype==0)=[];
    for ct = 1:length(ctype)
        if select_nochange==1
            selectepochs = intersect(no_change_conds,find(ep_mark==ctype(ct)));
        else
            selectepochs = find(ep_mark==ctype(ct));
        end
        EEG = pop_select(EEG_ALL,'trial',selectepochs);
        EEG.conditionlabel = strread(num2str(ep_mark(find(ep_mark==ctype(ct)))),'%s')'; 
        %sname = [C{1} '_cond' num2str(ctype(ct)) '.set'];
        %EEG = pop_saveset(EEG,'filename',sname,'filepath',raw_path); 
        if ALLEEG_save==1
            if ct==1
                ALLEEG = EEG; 
            else 
                ALLEEG(ct)=EEG;
            end
        end
    end
    if ALLEEG_save==1
        if select_nochange
            sname = [C{1} '_' C{2} '_nochange_ALLEEG.mat'];
        else
            sname = [C{1} '_' C{2} '_conds_ALLEEG.mat'];
        end
        save(sname,'ALLEEG');
    end
    
    if ALLEEG_save==2
        if f==files_ana(1)
            ALLEEG = EEG; 
        else 
            ALLEEG(length(ALLEEG)+1)=EEG;
        end
    end
end
if ALLEEG_save==2
    sname = ['Allfiles_nochange_ALLEEG.mat'];
    save(sname,'ALLEEG');
end