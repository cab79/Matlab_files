function varargout = SCIn(varargin)

% GUI:
% display messages stating trial number, block number, ready to startstop next block, etc.

%MATLAB code for SCIn.fig
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% Begin initialization code - DO NOT EDIT
dbstop if error

% set global variable d: list of directories
global d
d.root = 'C:\Matlab_files\NTIP\SCIn';
addpath(genpath(d.root))
d.expts = 'Experiments';
d.settings = 'Settings';
d.seq = 'Sequences';

gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @SCIn_OpeningFcn, ...
                   'gui_OutputFcn',  @SCIn_OutputFcn, ...
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


% --- Executes just before SCIn is made visible.
function SCIn_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to SCIn (see VARARGIN)

% Choose default command line output for SCIn
handles.output = hObject;
guidata(hObject, handles);

% set the window position
set(groot, 'Units', 'pixels') 
screensize = get(0, 'screensize');
set(hObject, 'Units', 'pixels')
set(hObject, 'Position', [screensize(3)/2 - 200,...
                          screensize(4)/2 - 200,...
                          400, 400]);

% --- Outputs from this function are returned to the command line.
function varargout = SCIn_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on selection change in ExptOpt.
function ExptOpt_Callback(hObject, eventdata, handles)
% hObject    handle to ExptOpt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns ExptOpt contents as cell array
%        contents{get(hObject,'Value')} returns selected item from ExptOpt
expts = get(hObject,'String');
opt = get(handles.ExptOpt,'Value');
[~,handles.exptFun,] = fileparts(expts{opt});
set(handles.info, 'String', 'Experiment selected');
guidata(hObject, handles)


% --- Executes during object creation, after setting all properties.
function ExptOpt_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ExptOpt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
global d
files = dir(fullfile(d.root,d.expts,'*.m'));
names = {files(:).name};
set(hObject,'String',names);

% --- Executes on selection change in SettingsOpt.
function SettingsOpt_Callback(hObject, eventdata, handles)
% hObject    handle to SettingsOpt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns SettingsOpt contents as cell array
%        contents{get(hObject,'Value')} returns selected item from SettingsOpt
settings = get(hObject,'String');
opt = get(handles.SettingsOpt,'Value');
[~,handles.SettingsFun,~] = fileparts(settings{opt});
set(handles.info, 'String', 'Settings selected');
guidata(hObject, handles)

% --- Executes during object creation, after setting all properties.
function SettingsOpt_CreateFcn(hObject, eventdata, handles)
% hObject    handle to SettingsOpt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
global d
files = dir(fullfile(d.root,d.settings,'*.m'));
names = {files(:).name};
set(hObject,'String',names);

% --- Executes on selection change in SeqOpt.
function SeqOpt_Callback(hObject, eventdata, handles)
% hObject    handle to SeqOpt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns SeqOpt contents as cell array
%        contents{get(hObject,'Value')} returns selected item from SeqOpt
seq = get(hObject,'String');
opt = get(handles.SeqOpt,'Value');
handles.SeqName = seq{opt};
set(handles.info, 'String', 'Sequence selected');
guidata(hObject, handles)

% --- Executes during object creation, after setting all properties.
function SeqOpt_CreateFcn(hObject, eventdata, handles)
% hObject    handle to SeqOpt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
global d
files = dir(fullfile(d.root,d.seq,'*.mat'));
names = {files(:).name};
set(hObject,'String',names);

% --- Executes on button press in CreateSeq.
function CreateSeq_Callback(hObject, eventdata, handles)
% hObject    handle to CreateSeq (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of CreateSeq

% load functions if not already loaded
set(handles.info, 'String', 'Loading settings...');
if ~isfield(handles,'SettingsFun')
    settings = get(handles.SettingsOpt,'String');
    opt = get(handles.SettingsOpt,'Value');
    [~,handles.SettingsFun,] = fileparts(settings{opt});
end
eval(['handles = ' handles.SettingsFun '(handles);']);
set(handles.info, 'String', 'Creating sequence...');
pause(0.1) % otherwise message not displayed
eval(['handles = ' handles.SeqFun '(handles);']);
set(handles.info, 'String', 'Sequence created. Save it?');
guidata(hObject, handles)

% --- Executes on button press in SaveSeq.
function SaveSeq_Callback(hObject, eventdata, handles)
% hObject    handle to SaveSeq (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hint: get(hObject,'Value') returns toggle state of SaveSeq
set(handles.info, 'String', 'Saving sequence...');
filename = ['AuditorySCIn_Dur_' num2str(handles.dur) '_f0_' num2str(handles.f0) '_df_' num2str(handles.df) '_deviant_' num2str(handles.fc)];
seq=handles.Seq;
global d
uisave({'seq'},fullfile(d.root,d.seq,filename));
set(handles.info, 'String', 'Sequence saved.');

% --- Executes on button press in StartStop.
function StartStop_Callback(hObject, eventdata, handles)
% hObject    handle to StartStop (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% check if the button is pressed
if get(hObject, 'Value') == get(hObject, 'Max')

    % info
    set(handles.info, 'String', 'Loading sequence...');

    % load functions if not already loaded
    if ~isfield(handles,'exptFun')
        expts = get(handles.ExptOpt,'String');
        opt = get(handles.ExptOpt,'Value');
        [~,handles.exptFun,] = fileparts(expts{opt});
    end
    if ~isfield(handles,'SeqName')
        seq = get(handles.SeqOpt,'String');
        opt = get(handles.SeqOpt,'Value');
        handles.SeqName = seq{opt};
    end
    if ~isfield(handles,'Seq')
        % load sequence
        global d
        load(fullfile(d.root,d.seq,handles.SeqName));
        handles.Seq = seq;
    end
    guidata(hObject, handles)


    % disable all GUI components except the toggle buttons
    set(handles.ExptOpt, 'Enable', 'off')
    set(handles.SeqOpt, 'Enable', 'off')
    set(handles.SettingsOpt, 'Enable', 'off')
    %set(handles.StartStop, 'Enable', 'off')

    % info
    set(handles.info, 'String', 'Running sequence...');
    
    % startstop
    opt = 'start';
    eval(['handles = ' handles.exptFun '(handles,opt);']);

    % change the button text
    set(hObject, 'String', 'Stop')
    
else

    % enable all GUI componets if the toggle buttons are unpressed
    if get(handles.StartStop, 'Value') == get(handles.StartStop, 'Min')
        set(handles.ExptOpt, 'Enable', 'on')
        set(handles.SeqOpt, 'Enable', 'on')
        set(handles.SettingsOpt, 'Enable', 'on')
        %set(handles.StartStop, 'Enable', 'on')
    else
        set(handles.ExptOpt, 'Enable', 'off')
        set(handles.SeqOpt, 'Enable', 'off')
        set(handles.SettingsOpt, 'Enable', 'off')
        %set(handles.StartStop, 'Enable', 'off')
    end

    % info
    set(handles.info, 'String', 'Stopped.');
    
    % stop 
    opt = 'stop';
    eval(['handles = ' handles.exptFun '(handles,opt);']);

    % change the button text
    set(hObject, 'String', 'Start')
end
    

% --- Executes on button press in PauseResume.
function PauseResume_Callback(hObject, eventdata, handles)
% hObject    handle to PauseResume (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of PauseResume

% check if the button is pressed
if get(hObject, 'Value') == get(hObject, 'Max')

    % info
    set(handles.info, 'String', 'Paused.');
    
    % pause
    opt = 'pause';
    eval(['handles = ' handles.exptFun '(handles,opt);']);

    % change the button text
    set(hObject, 'String', 'Resume')
    
else

    % info
    set(handles.info, 'String', 'Running sequence...');
    
    % resume 
    opt = 'resume';
    eval(['handles = ' handles.exptFun '(handles,opt);']);

    % change the button text
    set(hObject, 'String', 'Pause')
end
    
