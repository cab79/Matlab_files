function g=align_gplots(g,gap)

set(gcf, 'Units', 'normalized')
figpos=get(gcf,'Position');

for gi = 1:length(g)-1
    % remove xticks for all but last plot
    g(gi,1).facet_axes_handles.XTickLabel = [];
    g(gi,1).facet_axes_handles.XLabel = [];
end

g.draw();

% align ticks
txt='';
for gi = 1:length(g)
    txt = [txt 'g(' num2str(gi) ', 1).facet_axes_handles,'];
end
eval(['linkaxes([' txt '],''x'');']);

% align left and right
for gi = 1:length(g)
    pos(gi,:) = g(gi, 1).facet_axes_handles.Position;
    legpos(gi,:) = g(gi, 1).legend_axe_handle.Position;
    ylabpos(gi,:) = g(gi, 1).facet_axes_handles.YLabel.Position;
end
maxpos = max(pos,[],1);
minpos = min(pos,[],1);
minpos(3)=0.5;
minlegpos = min(legpos,[],1);
leg_start = (maxpos(1)+minpos(3));
leg_end = leg_start + minlegpos(3);
minylabpos = min(ylabpos,[],1);
figpos(3) = figpos(3)+figpos(3)*(((maxpos(1)-minpos(1))+(leg_end-(minlegpos(1)+minlegpos(3)))))*1.2;
set(gcf,'Units', 'normalized','Position',figpos);
for gi = 1:length(g)
    g(gi, 1).facet_axes_handles.Position(1) = maxpos(1);
    g(gi, 1).facet_axes_handles.Position(4) = maxpos(4)/gap; % remove gap size
    g(gi, 1).facet_axes_handles.Position(3) = minpos(3);
    g(gi, 1).legend_axe_handle.Position(1) = leg_start;
    g(gi, 1).facet_axes_handles.YLabel.Position(1) = minylabpos(1);
end

% set gap between subplots
for gi = length(g)-1:-1:1
    y_extent = g(gi+1, 1).facet_axes_handles.Position(4);
    g(gi, 1).facet_axes_handles.Position(2) = g(gi+1, 1).facet_axes_handles.Position(2)+y_extent*gap;
end

% set y position of legends
for gi = 1:length(g)
    y_extent = g(gi, 1).facet_axes_handles.Position(4);
    g(gi, 1).legend_axe_handle.Position(2) = g(gi, 1).facet_axes_handles.Position(2)+y_extent*0.05;
    g(gi, 1).legend_axe_handle.Title.Position(2) = g(gi, 1).legend_axe_handle.Title.Position(2)*0.05;
end
    