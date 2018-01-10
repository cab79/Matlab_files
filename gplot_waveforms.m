
%% run
function gplot_waveforms(P)
dbstop if error

%% load EEGLAB data to plot topography - will take a few mins if ERP_DAT file does not already exist
if P.plot_topo && ~exist(fullfile(P.eeglab_path,'ERP_DAT.mat'),'file')
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
    DAT=cell(length(fnames),length(P.eventtypes));
    for f = 1:length(fnames)
        EEGall=pop_loadset([fnames{f} '.set'],P.eeglab_path);
        for e = 1:length(P.eventtypes)
            if any(strcmp({EEGall.event.type},P.eventtypes{e}))
               EEG = pop_selectevent(EEGall,'type',P.eventtypes{e});
               DAT{f,e} = mean(EEG.data,3);
            end
        end
    end
    DAT=reshape(DAT,[],1);
    EEGtimes = EEGall.times;
    chanlocs=EEGall.chanlocs;
    save(fullfile(P.eeglab_path,'ERP_DAT.mat'),'DAT','chanlocs','EEGtimes');
elseif P.plot_topo
    load(fullfile(P.eeglab_path,'ERP_DAT.mat'));
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
if any(strcmp(P.fact,'Time')) && S.timelev
    time_ind = S.Fm(:,1+find(strcmp(P.fact,'Time')));
    fact_ind = fact_ind(time_ind==S.timelev,:);
    sub_ind = sub_ind(time_ind==S.timelev);
end

if ~isempty(P.plotclus)
    % restrict to selected clusters
    clnames = P.plotclus';
else
    % obtain all cluster names
    clnames = fieldnames(S.wf);
end

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
    % select wf if struct (e.g. source data)
    if isstruct(wf)
        itimes = wf.time*1000;
        wf=wf.wf;
    end
    % identify unique rows of the combination of factors
    [~,Frows,WFrows] = unique([fact_ind sub_ind],'rows');
    
    % create new factor indices for unique rows.
    Fi = fact_ind(Frows,:);
    Si = sub_ind(Frows,:);
    
    % average wf data over non-unique rows
    wff = cell(length(unique(WFrows)),1);
    for i = unique(WFrows)'
        if size(wf{1},1)>size(wf{1},2) % timepoints in first dimension
            wff{i}=mean(cell2mat(wf(WFrows==i)'),2);
        else
            wff{i}=mean(cell2mat(wf(WFrows==i)),1)';
        end
    end
    
    % if there is more than one factor, the first factor levels produce a
    % plot each
    ptitle={};
    Fi_ind={};
    if length(P.facplot)>1
        Nplots=length(P.cval{1,2});
        % titles
        for pt = 1:Nplots
            ptitle{pt} = P.cval{1,1}{P.cval{1,2}(pt)};
            Fi_ind{pt} = find(Fi(:,1)==P.cval{1,2}(pt)); % indices of Fi (and wf) for each plot
        end
    else
        Nplots = 1;
        ptitle{1} = '';
        Fi_ind{1} = 1:length(Fi);
    end
    
    % extract cluster statistics
    ct_ind = find(ismember(clustable(:,1),cllabel));
    if ~isempty(P.poly)
        E_val = P.poly;
    else
        %F_val = clustable{ct_ind,9};
        E_val = [clustable{ct_ind,11:12}];
        P_val = [clustable{ct_ind,10}];
    end
    
%end

% for each cluster, 
%for cl = 1:length(clnames)
    
    clear g
    for p = 1:Nplots
        % plot data trajectories
        % gramm supports 2D inputs for X and Y data (as 2D array or cell of
        % arrays), which is particularly useful for representing repeated
        % trajectories. The grouping data is then given per trajectory, given
        % as a 1xNtraj cellstr.
        
        % construct condition labels (of last or only factor)
        condlev=P.cval{end,1};
        factind=Fi(Fi_ind{p},end);
        fi{p}=ismember(factind,P.cval{end,2});
        factind = factind(fi{p});
        cond{p} = condlev(factind);
        
        % set x and y axis data to plot
        y=wff(Fi_ind{p});
        y=y(fi{p});
        if ~isempty(P.xlimits)
            xlim = dsearchn(itimes',P.xlimits')';% x values of the selected segment to plot
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
        g(1,p)=gramm('x',x,'y',y,'color',cond{p});

        % name the axes and legend
        g(1,p).set_names('x',P.xaxisname,'y',P.yaxisname,'color',P.fact_names{end},'column','','row','');

        g(1,p).set_order_options('color',-1);
        g(1,p).set_color_options('map',P.colours);
        g(1,p).set_point_options('base_size',3);
        g(1,p).set_title(ptitle{p});
        g(1,p).geom_vline('xintercept',P.xlines,'style','k--','extent',4);
        g(1,p).geom_vline('xintercept',P_val,'style','k');
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
    
    if P.save_waveforms
        psize = 18;
        g.export('file_name',[cllabel '_WF'],'export_path',fullfile(P.spm_path,P.clusdir),'file_type','png','width',Nplots*psize,'height',psize,'units','centimeters');
    end
    
    if P.plot_topo
        % average ERPs over non-unique rows
        DATa = cell(length(unique(WFrows)),1);
        for i = unique(WFrows)'
            DATa{i}=mean(cat(3,DAT{WFrows==i}),3);
        end
        f1=figure
        pln=0;
        for p = 1:Nplots
            y=DATa(Fi_ind{p});
            y=y(fi{p});
            [~,~,condind]=unique(cond{p},'stable');
            unicond = unique(condind)';
            for cn = unicond
                peakdata=y(condind==cn);
                peakdata = mean(cat(3,peakdata{:}),3);
                plotchans=1:length(chanlocs);
                plotchans(P.no_plot_ele)=[];
                %[~,markchans] = intersect(plotchans,tp);
                if any(P.topo_subtimewin) && length(P.topo_subtimewin)==1 % multiple plots within a range
                    cluswin = E_val(end)-E_val(1);
                    nlat = floor(cluswin/P.topo_subtimewin);
                    if nlat>1
                        lats = linspace(E_val(1),E_val(end),nlat);
                        lats = dsearchn(EEGtimes',lats');
                    else
                        lats = [find(EEGtimes==E_val(1)) find(EEGtimes==E_val(end))];
                    end
                    
                    for ln = 1:length(lats)-1
                        pln = pln+1;
                        lat=lats(ln:ln+1);
                        subplot(Nplots+length(unicond),length(lats)-1,pln); topoplot(mean(peakdata(:,lat(1):lat(2)),2), chanlocs,'maplimits','absmax','electrodes','on','plotchans',plotchans);%,'emarker2',{markchans,'o','w',7,1}); 
                        title(num2str(EEGtimes(lat)'))
                    end
                elseif any(P.topo_subtimewin) && length(P.topo_subtimewin)==2 % specified time window
                    pln = pln+1;
                    lat = dsearchn(EEGtimes',P.topo_subtimewin');
                    subplot(Nplots+length(unicond),Nplots,pln); topoplot(mean(peakdata(:,lat),2), chanlocs,'maplimits','absmax','electrodes','on','plotchans',plotchans);%,'emarker2',{markchans,'o','w',7,1}); 
                    title(num2str(EEGtimes(lat)'))
                else
                    pln = pln+1;
                    lat = find(EEGtimes==E_val(1));
                    subplot(Nplots+length(unicond),Nplots,pln); topoplot(mean(peakdata(:,lat),2), chanlocs,'maplimits','absmax','electrodes','on','plotchans',plotchans);%,'emarker2',{markchans,'o','w',7,1}); 
                    title(num2str(EEGtimes(lat)'))
                end
            end
            %colorbar
        end
        tightfig(f1)
        %save figure 
        if P.save_topo
            print(fullfile(P.spm_path,P.clusdir,[cllabel '_topo']),'-dpng');
        end
    end
        
end
%x=range
end