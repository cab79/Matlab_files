function nested_callback_ex
% Create and configure timer object
t = timer('ExecutionMode','fixedRate', ... % Run continuously
'TimerFcn',@MyTimerFcn); % Run MyTimerFcn at each timer event
% Load and display a grayscale image
Im = imread('street1.jpg');imagesc(Im);
%i=0;
% Start the timer
start(t)
pause(2)
stop(t)
pause(2)
start(t)
%get(t, 'AveragePeriod')


function MyTimerFcn(obj,event)
% Scale and display image
i=1;
while strcmp(get(t,'Running'),'on')
    Im = Im*1.1; % Make brighter
    imagesc(Im) % Display updated image
    %get(t, 'InstantPeriod')
    i = (i-1)+1
    pause(1)
end
end


end % function nested_callback_ex



