function g=align_gplots(g,gap)

% remove yticks for all but last plot
for gi = 1:length(g)-1
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
minlegpos = min(legpos,[],1);
minylabpos = min(ylabpos,[],1);
for gi = 1:length(g)
    g(gi, 1).facet_axes_handles.Position(1) = maxpos(1)*1.2;
    g(gi, 1).facet_axes_handles.Position(4) = maxpos(4)/gap; % remove gap size
    g(gi, 1).facet_axes_handles.Position(3) = minpos(3);
    g(gi, 1).legend_axe_handle.Position(1) = minlegpos(1) + maxpos(1)*0.2;
    g(gi, 1).facet_axes_handles.YLabel.Position(1) = minylabpos(1);
end

% set gap between subplots
for gi = length(g)-1:-1:1
    y_extent = g(gi+1, 1).facet_axes_handles.Position(4);
    g(gi, 1).facet_axes_handles.Position(2) = g(gi+1, 1).facet_axes_handles.Position(2)+y_extent*gap;
end

% set position of legends
for gi = 1:length(g)
    y_extent = g(gi, 1).facet_axes_handles.Position(4);
    g(gi, 1).legend_axe_handle.Position(2) = g(gi, 1).facet_axes_handles.Position(2)+y_extent*0.05;
    g(gi, 1).legend_axe_handle.Title.Position(2) = g(gi, 1).legend_axe_handle.Title.Position(2)*0.05;
end
    