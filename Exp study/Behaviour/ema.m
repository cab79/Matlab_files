function out = ema(data,period)
% Function to calculate the exponential moving average of a data set
% 'data' is the vector to operate on.  The first element is assumed to be
% the oldest data.
% 'period' is the number of periods over which to calculate the average
%
% Example:
% out = ema(data,period)
%

% Error check
if nargin ~= 2
    error([mfilename,' requires 2 inputs.']);
end
[m,n]=size(data);
if ~(m==1 || n==1)
    error(['The data input to ',mfilename,' must be a vector.']);
end
if (numel(period) ~= 1)
    error('The period must be a scalar.');
end

% convert the period to an exponential percentage
ep = 2/(period+1);
% calculate the EMA
out = filter(ep,[1 -(1-ep)],data,data(1)*(1-ep));
out(1:period-1)=nan;