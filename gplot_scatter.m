function g = gplot_scatter(g,P)

%construct gramm plot
g(P.xy(1),P.xy(2))=gramm('x',P.x,'y',P.y,'color',P.cond,'size',P.condsize);
g(P.xy(1),P.xy(2)).set_names('x',P.xaxisname,'y',P.yaxisname,'color','','column','','row','');
%g(P.xy(1),P.xy(2)).set_order_options('color',-1);
g(P.xy(1),P.xy(2)).set_color_options('map',P.colours);
g(P.xy(1),P.xy(2)).set_point_options('base_size',1);
g(P.xy(1),P.xy(2)).set_line_options('base_size',2);
g(P.xy(1),P.xy(2)).set_title(P.ptitle);
if ~P.legend
    g(P.xy(1),P.xy(2)).no_legend();
end

switch P.plottype 
    case 'geom_point'
        g(P.xy(1),P.xy(2)).geom_point();
    case 'geom_point_line'
        g(P.xy(1),P.xy(2)).geom_point();
        g(P.xy(1),P.xy(2)).stat_glm('geom','line');

end
end
%tix=get(g(P.xy(1),P.xy(2)).facet_axes_handles,'ytick')
%g(P.xy(1),P.xy(2)).axe_property('yticklabel',num2str(tix,'%.1f'));

%Possibility to set color and fill by indices (using a column vector of
%integers. Colormap generated between 1 and max(vector))
%   g(2,1).geom_polygon('y',{[5 20];  [20 30];  [30 50]},'color',[1 ; 3;  2]);
