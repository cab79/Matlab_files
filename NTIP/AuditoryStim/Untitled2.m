function test
t=timer;
set(t,'ExecutionMode','fixedDelay','BusyMode','drop','Period',0.001); 
t.TimerFcn = @my_callback_fcn; 
start(t); 
stop(t); 
get(t, 'AveragePeriod')

function my_callback_fcn(obj,event)
    x=1;
end

end