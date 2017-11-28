% Copyright 2010 - The MathWorks, Inc.
%% Using Data Acquisition Toolbox to write 24 lines of digital output data
% to a National Instruments NI-USB-6501 device

dio=digitalio('nidaq', 'Dev1'); % create digital IO object
addline(dio, 0:23, 'out'); % add 24 lines representing the 3 ports of 8 lines each
lineconfig=dio.Line % display index and digital line configurations to show mapping

%% Output the data using binvec (binary vector)

% Note: the code below will change the digital output voltage of your DIO device 
putvalue(dio.Line(1:8), [0 0 0 1 1 0 0 1]);%Port 0.0 through 0.7 
putvalue(dio.Line(9:16), [1 1 1 1 1 1 1 1]);%Port 1.0 through 1.7
putvalue(dio.Line(17:24), [1 1 1 1 0 0 0 0]);%Port 2.0 through 2.7

%% Convert binvec to decimal

% MATLAB converts the binvec to a decimal value(0-255); value is displayed
% for illustration 
% binvec2dec is a Data Acquisition Toolbox command
byte1=binvec2dec([0 0 0 1 1 0 0 1])% byte1=152
byte2=binvec2dec([1 1 1 1 1 1 1 1])% byte2=255
byte3=binvec2dec([1 1 1 1 0 0 0 0])% byte3=15

%% Output the data using decimal values (generates same output as putvalue statements above)

putvalue(dio.Line(1:8),byte1); 
putvalue(dio.Line(9:16),byte2); 
putvalue(dio.Line(17:24),byte3);



