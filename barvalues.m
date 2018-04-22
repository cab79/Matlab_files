function varargout = barvalues(h,precision)
% barvalues;
% barvalues(h)
% barvalues(h,precision)
% t = barvalues(h,precision)
%
% Display bar values ontop of bars in bar plot.
% 
% h - handle to axes or bar (operates on specified object only) 
%     or figure (operates on all desndant bars).
%     (default:gca)
% precision - Decimal precision to display (0-10),
%             or 'formatSpec' as in num2str. (default:'% .0f')
% t - handles to the text objects.

%Author: Elimelech Schreiber, 11/2017 
% ver 1.3

t=[];

if nargin>1 && ~isempty(precision) % Parse precision
    if isnumeric(precision) && precision >=0 && precision <=10
        precision =['% .',int2str(precision),'f'];
    elseif ~ischar(precision) && ~isstring(precision)
        error('Precision format unsupported.');
    end
else
    precision ='% .0f';
end

if nargin<1 || isempty(h)   % parse h (handle)
    h =gca;
elseif isaType(h,'figure')
   B =findobj(h,'type','bar'); % apply to multiple axes in figure.
   for b =B'
           t = [t; {barvalues(b,precision)}]; % Return array of text objects
                                              % for each bar plot.
   end
    if nargout>0
        varargout{1}=t;
    end
    return;
end
if isaType(h,'axes')
    h =findobj(h,'type','bar');
    if isempty(h)
        return; % silently. to support multiple axes in figure.
    end
end
if ~isaType(h,'bar')
    error('Cannot find bar plot.');
end

for hn =h'
    
    axes(ancestor(hn,'axes')); % make intended axes curent.
    if isfield(hn,'XOffset')&&~isempty(hn.XOffset), XOffset = hn.XOffset; else XOffset = 0; end
    if isfield(hn,'YOffset')&&~isempty(hn.YOffset), YOffset = hn.YOffset; else YOffset = 0; end
    xData = hn.XData +XOffset; yData = hn.YData +YOffset;
     
    t = [t;  text(xData,yData,...    %position
       arrayfun(@(x)num2str(x,precision),yData,'UniformOutput' ,false),...    %text to display
        'HorizontalAlignment','center','VerticalAlignment','bottom')];
end
if nargout>0
    varargout{1}=t;
end

function flag =isaType(h,type)
try
    flag =strcmp(get(h, 'type'), type); 
catch
    flag =false;
end


function flag = isfield(h,fld)
flag =true;
try
    get(h,fld)
catch
    flag =false;
end
