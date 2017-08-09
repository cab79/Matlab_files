function plotdotbox(P) 

nf = regexp(P.filenames,'*');
if ~isempty(nf{:})
    files = dir(fullfile(P.pth,P.filenames{:}));
    fnames = {files(:).name};
    datnames=repmat(P.datnames,1,length(fnames));
    yaxisnames=repmat(P.yaxisnames,1,length(fnames));
    xaxisnames =repmat(P.xaxisnames,1,length(fnames));
    subplotnames=repmat(P.subplotnames,1,length(fnames));
    groupingnames =repmat(P.groupingnames,1,length(fnames));
else
    fnames = P.filenames;
    datnames=P.datnames;
    yaxisnames=P.yaxisnames;
    xaxisnames =P.xaxisnames;
    subplotnames=P.subplotnames;
    groupingnames =P.groupingnames;
end

% plot as separate figures for each file?
if P.plot_sep_fg == 1
    Npt = length(fnames);
    Nfg = 1;
    fignames = {''};
    pltnames = fnames;
else
    Npt = 1;
    Nfg = length(fnames);
    fignames = fnames;
    pltnames = {''};
end

plottypes = P.plottypes([P.plottypes{:,2}]==1,1); 

for fg = 1:Nfg
    clear g
    for f = 1:Npt
        fi = max(f,fg);
        
        dat = xls2struct(fullfile(P.pth,fnames{fi}));
        %dat.Data = log(dat.Data+sqrt(dat.Data.^2+1)); %arsinh
        %dat.Data = sign(dat.Data).*log(abs(dat.Data)+1); %log-modulus
        
        % extract data from each cell for normality tests
        %[datnum,dattxt,datraw] = xlsread(fullfile(P.pth,fnames{fi}));
        %factdat=dattxt(2:end,1:end-2);
        %for fd = 1:size(factdat,2)
        %    [~,~,ind] = unique(factdat(:,fd));
        %    factdatnum(:,fd) = ind;
        %end
        %[~,~,cells] = unique(factdatnum,'rows','stable');
        %for c=1:length(unique(cells))
        %    %h(c,fi) = kstest(dat.Data(cells==c));
        %    [h(c,fi),pval(c,fi),~] = swtest(dat.Data(cells==c), 0.05);
        %end
        
        
        % Change data format for line plots
        if any(strcmp(plottypes,'line'))
            xvar = unique(dat.(genvarname(xaxisnames{fi})),'stable');
            xvar = 1:length(xvar);
            [Sn,Su,Si]=unique(dat.Subject,'stable');
            ydata = dat.(genvarname(datnames{fi}));
            if ~isempty(subplotnames)
                sublevels = dat.(genvarname(subplotnames{fi}));
                if isnumeric(P.sublevel)
                    yi = find(sublevels == P.sublevel);
                end
                Si=Si(yi);
                ydata = ydata(yi);
            end
            for s = 1:length(Sn)
                yvar{s,1}=ydata(Si==s);
            end
            gvar = dat.(genvarname(groupingnames{fi}));
            gvar = gvar(Su)';
            lvar = Sn;
            g(f,1)=gramm('x',xvar,'y',yvar,'color',gvar);%,'marker',dat.Group);%,'subset',cars.Cylinders~=3 & cars.Cylinders~=5);
        else
            xvar = dat.(genvarname(xaxisnames{fi}));
            yvar = dat.(genvarname(datnames{fi}));
            gvar = dat.(genvarname(groupingnames{fi}));
            lvar = dat.Subject;
            g(f,1)=gramm('x',xvar,'y',yvar,'color',gvar,'label',lvar);%,'marker',dat.Group);%,'subset',cars.Cylinders~=3 & cars.Cylinders~=5);
            % Subdivide the data in subplots
            if ~isempty(subplotnames{fi})
                g(f,1).facet_grid([],dat.(genvarname(subplotnames{fi})));
            end
        end
        

        % names the axes and legend
        g(f,1).set_names('x',xaxisnames{fi},'y',yaxisnames{fi},'color',groupingnames{fi},'column','','row','');
        g(f,1).set_title(pltnames{f});

        template = g(f,1);
        
        p = 0;
        %Jittered scatter plot
        if any(strcmp(plottypes,'jitter'))
            p=p+1;
            g(f,p)=copy(template);
            g(f,p).geom_jitter('width',0.4,'height',0);
            if P.label_points; g(f,1).geom_label('color','k','dodge',0.7,'VerticalAlignment','bottom','HorizontalAlignment','center');end
        end
        %g(f,1).set_title('geom_jitter()');

        %Boxplots
        if any(strcmp(plottypes,'boxplot'))
            p=p+1;
            g(f,p)=copy(template);
            g(f,p).stat_boxplot();
        end
        %g(f,2).set_title('stat_boxplot()');
        
        %Violin plots
        if any(strcmp(plottypes,'violin'))
            p=p+1;
            g(f,p)=copy(template);
            g(f,p).stat_violin('fill','transparent');
        end

        %Averages with confidence interval
        if any(strcmp(plottypes,'confidence'))
            p=p+1;
            g(f,p)=copy(template);
            g(f,p).stat_summary('geom',{'bar','black_errorbar'});
        end
        
        %Line plots
        if any(strcmp(plottypes,'line'))
            p=p+1;
            g(f,p)=copy(template);
            g(f,p).geom_line('dodge',0,'alpha',1);
        end

        %Raw data as scatter plot
        %g(2,1).geom_point();
        %g(2,1).set_title('geom_point()');

        %These functions can be called on arrays of gramm objects
        %g.set_title('Visualization of Y~X relationships with X as categorical variable');
    end
    g.set_title(fignames{fg});
    g.set_text_options('base_size',P.textsize);
    if ~isempty(P.groupcolours)
        g.set_color_options('map',P.groupcolours);
    end
    %g.set_point_options('markers',{'o','s'});
    %figure('Position',[100 100 800 550]);
    fig=figure;%('Position',[100 100 800 550]);
    set(fig, 'Units', 'normalized', 'Position', [0,0,1,1]);
    g.draw();
    if P.save_figure
        sizefig = 18;
        C=strsplit(fnames{fi},'.');
        g.export('file_name',[C{1} '_plots'],'export_path',P.pth,'file_type','png','width',p*sizefig,'height',sizefig,'units','centimeters');
    end
end
%norm_head = horzcat({''},fnames);
%norm_rows = num2cell(unique(cells));
%h = vertcat(norm_head,horzcat(norm_rows,num2cell(h)))
%pval = vertcat(norm_head,horzcat(norm_rows,num2cell(pval)))