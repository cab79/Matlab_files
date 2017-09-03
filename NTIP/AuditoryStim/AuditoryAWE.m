function varargout = AuditoryAWE(varargin)
%MATLAB code for AuditoryAWE.fig
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% Begin initialization code - DO NOT EDIT
dbstop if error
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @AuditoryAWE_OpeningFcn, ...
                   'gui_OutputFcn',  @AuditoryAWE_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before AuditoryAWE is made visible.
function AuditoryAWE_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to AuditoryAWE (see VARARGIN)

% Choose default command line output for AuditoryAWE
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% set the window position
set(groot, 'Units', 'pixels') 
screensize = get(0, 'screensize');
set(hObject, 'Units', 'pixels')
set(hObject, 'Position', [screensize(3)/2 - 200,...
                          screensize(4)/2 - 200,...
                          400, 400]);
                      

% --- Outputs from this function are returned to the command line.
function varargout = AuditoryAWE_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


function finput_Callback(hObject, eventdata, handles)
% hObject    handle to finput (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% error check
f0 = str2double(get(hObject,'String'));
if isreal(f0) && ~isnan(f0) && ~isempty(f0) && f0 >= 0 && f0 <= 20000
    % do nothing
else
    errordlg('The tone frequency must be a real number in the interval 0 - 20000!',...
             'Tone Frequency Error');
    set(hObject, 'String', '300')
    return
end

% declare the variables f0 and fs in the Handles Structure
handles.fs = 96000;
handles.f0 = f0;
guidata(hObject, handles)


% --- Executes during object creation, after setting all properties.
function finput_CreateFcn(hObject, eventdata, handles)
% hObject    handle to finput (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function dfinput_Callback(hObject, eventdata, handles)
% hObject    handle to dfinput (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of dfinput as text
%        str2double(get(hObject,'String')) returns contents of dfinput as a double

% error check
df = str2double(get(hObject,'String'));
if isreal(df) && ~isnan(df) && ~isempty(df) && df >= 0 && df <= 20000
    % do nothing
else
    errordlg('The difference frequency must be a real number in the interval 0 - 20000!',...
             'Difference Frequency Error');
    set(hObject, 'String', '10')
    return
end

% declare the variable df in the Handles Structure
handles.df = df;
guidata(hObject, handles)


% --- Executes during object creation, after setting all properties.
function dfinput_CreateFcn(hObject, eventdata, handles)
% hObject    handle to dfinput (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in Start.
function Start_Callback(hObject, eventdata, handles)
% hObject    handle to Start (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

persistent LRchan
   
% collect the user defined parameters via
% execution of the corresponding callbacks
finput_Callback(handles.finput, eventdata, handles)
handles = guidata(hObject);
dfinput_Callback(handles.dfinput, eventdata, handles)
handles = guidata(hObject);
Dur_Callback(handles.Dur, eventdata, handles)
handles = guidata(hObject);

% check if the button is pressed
if get(hObject, 'Value') == get(hObject, 'Max')
    
    t = transpose((1:handles.dur*handles.fs)/handles.fs);
    
    % construct the player object: left
    x = sin(2*pi*handles.f0*t);
    % construct the player object: right
    y = sin(2*pi*(handles.f0+handles.df)*t);
    
    % ERP / freq tag settings
    fpitch = 0; % Hz - 1/fpitch must be an integer multiple of 1/handles.df 
    pitchdiff = 100;
    finten = 2.5; % Hz - 1/finten must be an integer multiple of 1/handles.df 
    intendiff = 0.75; % multiple of normal intensity (value between 0 and 1)
    
    % pitch changes
    if fpitch>0
        % alternate sin waves
        x2 = sin(2*pi*(handles.f0+pitchdiff)*t);
        y2 = sin(2*pi*(handles.f0+handles.df+pitchdiff)*t);
        % define durations
        i1 = ones((1/fpitch)*handles.fs,1); 
        i2 = zeros((1/fpitch)*handles.fs,1); 
        i12 = repmat([i1;i2],handles.dur*(handles.fs/(length(i1)+length(i2))),1);
        % splice them in
        x(find(i12)) = x2(find(i12));
        y(find(i12)) = y2(find(i12));
    end
    
    % intensity changes
    if finten>0
        % alternate sin waves
        x2 = intendiff*sin(2*pi*(handles.f0)*t);
        y2 = intendiff*sin(2*pi*(handles.f0+handles.df)*t);
        % define durations
        i1 = ones((1/finten)*handles.fs,1); 
        i2 = zeros((1/finten)*handles.fs,1); 
        i12 = repmat([i1;i2],handles.dur*(handles.fs/(length(i1)+length(i2))),1);
        % splice them in
        x(find(i12)) = x2(find(i12));
        y(find(i12)) = y2(find(i12));
    end
    
    LRchan = audioplayer([x y], handles.fs);
    
    % disable all GUI componets except the toggle buttons
    set(handles.finput, 'Enable', 'off')
    set(handles.dfinput, 'Enable', 'off')
    set(handles.Dur, 'Enable', 'off')
    
    % play the sine-wave tone on the left audio channel
    play(LRchan)
    
    % change the button text
    set(hObject, 'String', 'Stop')
    
else
    
    % enable all GUI componets if the toggle buttons are unpressed
    if get(handles.Start, 'Value') == get(handles.Start, 'Min')
        set(handles.finput, 'Enable', 'on')
        set(handles.dfinput, 'Enable', 'on')
        set(handles.Dur, 'Enable', 'on')
    else
        set(handles.finput, 'Enable', 'off')
        set(handles.dfinput, 'Enable', 'off')
        set(handles.Dur, 'Enable', 'off')
    end
    
    % stop the left channel
    stop(LRchan)
    
    % change the button text
    set(hObject, 'String', 'Start')
    
end


% --- Executes on slider movement.
function slider1_Callback(hObject, eventdata, handles)
% hObject    handle to slider1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider


% --- Executes during object creation, after setting all properties.
function slider1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end




function Dur_Callback(hObject, eventdata, handles)
% hObject    handle to Dur (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Dur as text
%        str2double(get(hObject,'String')) returns contents of Dur as a double

Dur = str2double(get(hObject,'String'));
handles.dur = Dur;
guidata(hObject, handles)

% --- Executes during object creation, after setting all properties.
function Dur_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Dur (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
