%% plot ERP waveforms from clusters using the amazing Gramm toolbox!

% PREREQUISITES: 
% - An SPM.mat file containing a contrast on sensor-space ERP data 
% - Associated cluster data created from "Extract_clusters.m"
% - Cluster data generated from "Extract_cluster_waveforms.m"

clear all
close all
%% generic directories for all analyses for this study
%-------------------------------------------------------------
% directory in which SPM is saved
P.spm_path = 'C:\Data\Catastrophising study\SPMstats\Include1\Between\t-5500_-2500_b-5500_-5000_m_Grp_Att_Exp_Subject_orig_cleaned_spm';

%cluster directory name, which also specifies the constrast that will be
%plotted (i.e. the characters before the underscore)
P.clusdir='Grp_Exp_clusters';

%% specific directory and file information for this analysis
%-------------------------------------------------------------
%generic cluster waveform file name
P.wfname = 'cluster_data.mat';
%name of batch .mat file saved from design_batch.m and within same folder
%as SPM.mat
P.batch = 'matlabbatch.mat';
%name of 'subject' factor in the SPM design
P.subfactname = 'Subject';

%% plot options
%-------------------------------------------------------------
%CURRENTLY DOES NOT SUPPORT PLOTTING MORE THAN TWO FACTORS

%factor(s) to plot - if more than one, first factor levels will be on separate plots
%must be the same characters as used to name the factors in design_batch
%P.facplot={'Att','Exp'};
P.facplot={'Exp','Grp'};
%Full factor names to use for labelling the plots
fact_names = {
    %'Attention task';
    'Expectancy condition'
    'Group';
    };
%condition labels, i.e. levels of each condition, in the same order as in
%the SPM design matrix. One row per factor. second column in plotting order
cval={
    {'Cue 1 Low, Cue 2 Low','Cue 1 High, Cue 2 Low'};
    {'Low PCS','High PCS'}
    %{'No Task','Task'};
    };
colours = [0.2 0.5 1; 1 0.2 0.2]; % blue, red
% axis names - can differ according the experiment
P.xaxisname = {'peri-stimulus time (ms)'};
P.yaxisname = {'amplitude (arbitrary units)'};
% vertical dashed lines to indicate events, time in ms
xlimits = [];
% vertical dashed lines to indicate events, time in ms
xlines = [-5000 -2500 0];

% plot topographies
plot_topo=0;
eventtypes = {'c1','c2','c3','c4','c5','c6','c7','c8'};
eeglab_path = 'C:\Data\Catastrophising study\Preprocessed';
no_plot_ele=[];

% save optons
save_waveforms=0;
save_topo=0;


%% run

%% load EEGLAB data to plot topography - will take a few mins if ERP_DAT file does not already exist
if plot_topo && ~exist(fullfile(eeglab_path,'ERP_DAT.mat'),'file')
    % load SPM
    load(fullfile(P.spm_path,'SPM.mat'));
    imglist = SPM.xY.P; % Subject-condition image list
    st_string='mspm12_';
    en_string='\scond';
    fnames={};
    for i = 1:length(imglist);
        st=strfind(imglist{i},st_string);
        en=strfind(imglist{i},en_string);
        fnames{i,1} =imglist{i}(st+length(st_string):en-1);
    end
    [fnames,~,fni] = unique(fnames,'stable');
    DAT=cell(length(fnames),length(eventtypes));
    for f = 1:length(fnames)
        EEGall=pop_loadset([fnames{f} '.set'],eeglab_path);
        for e = 1:length(eventtypes)
            if any(strcmp({EEGall.event.type},eventtypes{e}))
               EEG = pop_selectevent(EEGall,'type',eventtypes{e});
               DAT{f,e} = mean(EEG.data,3);
            end
        end
    end
    DAT=reshape(DAT,[],1);
    EEGtimes = EEGall.times;
    save(fullfile(eeglab_path,'ERP_DAT.mat'),'DAT','chanlocs','EEGtimes');
elseif plot_topo
    load(fullfile(eeglab_path,'ERP_DAT.mat'));
end

%% plot cluster waveforms
% load waveform data
load(fullfile(P.spm_path,P.clusdir,P.wfname));

% load SPM design and list factors
load(fullfile(P.spm_path,P.batch));
P.fact = {matlabbatch{1,1}.spm.stats.factorial_design.des.fblock.fac(:).name};

% identify file indices relating to factor of interest and subject
fact_col = [];
for c = 1:length(P.facplot) % using a loop ensures that fact_col is in correct order for condition labels applied later
    fact_col(c) = find(ismember(P.fact,P.facplot(c)));
end
sub_col = find(ismember(P.fact,P.subfactname));
fact_ind = S.Fm(:,1+fact_col);
sub_ind = S.Fm(:,1+sub_col);

% obtain cluster names
clnames = fieldnames(S.wf);

% load example cluster image to extract time information
Cnii = load_nii(fullfile(P.spm_path,P.clusdir,[clnames{1} '.nii']));
t_res = Cnii.hdr.hist.srow_z(3);
t_off = Cnii.hdr.hist.srow_z(4);
t_dim = Cnii.hdr.dime.dim(4);
t_last = t_off+t_res*(t_dim-1);
itimes = t_off:t_res:t_last;

% load cluster statistics
load(fullfile(P.spm_path,P.clusdir,'cluster_table.mat'));

% colormap of background polygons
%colormap(gray);
%cmap=colormap;
%close all
range=[];
% for each cluster, 
for cl = 1:length(clnames)
    C=strsplit(clnames{cl},'_');
    cllabel = C{1};
    
    % select cluster data
    wf = S.wf.(clnames{cl});
    % identify unique rows of the combination of factors
    [~,Frows,WFrows] = unique([fact_ind sub_ind],'rows');
    
    % create new factor indices for unique rows.
    Fi = fact_ind(Frows,:);
    Si = sub_ind(Frows,:);
    
    % average wf data over non-unique rows
    wff = cell(length(unique(WFrows)),1);
    for i = unique(WFrows)'
        wff{i}=mean(cell2mat(wf(WFrows==i)'),2);
    end
    
    % if there is more than one factor, the first factor levels produce a
    % plot each
    ptitle={};
    Fi_ind={};
    if length(P.facplot)>1
        Nplots=max(fact_ind(:,1));
        % titles
        for pt = 1:Nplots
            ptitle{pt} = cval{1,1}{pt};
            Fi_ind{pt} = find(Fi(:,1)==pt); % indices of Fi (and wf) for each plot
        end
    else
        Nplots = 1;
        ptitle{1} = '';
        Fi_ind{1} = 1:length(Fi);
    end
    
    % extract cluster statistics
    ct_ind = find(ismember(clustable(:,1),cllabel));
    F_val = clustable{ct_ind,9};
    E_val = clustable{ct_ind,8};
    
    clear g
    for p = 1:Nplots
        % plot data trajectories
        % gramm supports 2D inputs for X and Y data (as 2D array or cell of
        % arrays), which is particularly useful for representing repeated
        % trajectories. The grouping data is then given per trajectory, given
        % as a 1xNtraj cellstr.
        
        % construct condition labels (of last or only factor)
        cond = cval{end,1}(Fi(Fi_ind{p},end));
        
        % set x and y axis data to plot
        y=wff(Fi_ind{p});
        if ~isempty(xlimits)
            xlim = dsearchn(itimes',xlimits')';% x values of the selected segment to plot
            x = itimes(xlim(1):xlim(2));
            for n = 1:length(y)
                y{n} = y{n}(xlim(1):xlim(2));
            end
        else
            x=itimes;
        end
        
        %change order (plot last

        %construct gramm plot
        %cond = cond(end:-1:1);
        %y = y(end:-1:1);
        g(1,p)=gramm('x',x,'y',y,'color',cond);

        % name the axes and legend
        g(1,p).set_names('x',P.xaxisname,'y',P.yaxisname,'color',fact_names{end},'column','','row','');

        g(1,p).set_order_options('color',-1);
        g(1,p).set_color_options('map',colours);
        g(1,p).set_point_options('base_size',3);
        g(1,p).set_title(ptitle{p});
        g(1,p).geom_vline('xintercept',xlines,'style','k--','extent',4);
        g(1,p).geom_vline('xintercept',E_val(1),'style','k');
        g(1,p).geom_polygon('x',{[min(E_val) max(E_val)]},'color',[0.5 0.5 0.5]);
        range(cl,:) = [min(E_val) max(E_val)];
        
        g(2,p) = copy(g(1,p));
        
        g(1,p).stat_smooth();
        g(1,p).set_line_options('base_size',0.5);
        g(2,p).stat_summary('type','ci','geom','area','setylim',true);
        

        %Possibility to set color and fill by indices (using a column vector of
        %integers. Colormap generated between 1 and max(vector))
     %   g(2,1).geom_polygon('y',{[5 20];  [20 30];  [30 50]},'color',[1 ; 3;  2]);

    end
    g.set_title(clnames{cl});
    g.set_text_options('base_size',12);
    fig=figure;%('Position',[100 100 800 550]);
    set(fig, 'Units', 'normalized', 'Position', [0,0,1,1]);
    g.draw();
    
    if save_waveforms
        size = 18;
        g.export('file_name',[cllabel '_WF'],'export_path',fullfile(P.spm_path,P.clusdir),'file_type','png','width',Nplots*size,'height',size,'units','centimeters');
    end
    
    if plot_topo
        % average ERPs over non-unique rows
        DATa = cell(length(unique(WFrows)),1);
        for i = unique(WFrows)'
            DATa{i}=mean(cat(3,DAT{WFrows==i}),3);
        end
        figure
        for p = 1:Nplots
            peakdata = mean(cat(3,DATa{Fi_ind{p}}),3);
            plotchans=1:length(chanlocs);
            plotchans(no_plot_ele)=[];
            %[~,markchans] = intersect(plotchans,tp);
            lat = find(EEGtimes==E_val(1));
            subplot(1,Nplots,p); topoplot(mean(peakdata(:,lat),2), chanlocs,'maplimits','absmax','electrodes','on','plotchans',plotchans);%,'emarker2',{markchans,'o','w',7,1}); 
            colorbar
        end
        %save figure 
        if save_topo
            print(fullfile(P.spm_path,P.clusdir,[cllabel '_topo']),'-dpng');
        end
    end
        
end
range