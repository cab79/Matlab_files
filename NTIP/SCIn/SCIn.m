function varargout = SCIn(varargin)

% GUI:
% display messages stating trial number, block number, ready to startstop next block, etc.
% select blocks with GUI by adding block number to h.Seq, e.g. row 4, and
% remove code from TSOT

%MATLAB code for SCIn.fig
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% Begin initialization code - DO NOT EDIT
%dbstop if error


global d
try
    rootdir
catch
    %error('start SCIn from the SCIn directory')
    error('update "rootdir.m" with the SCIn directory')
end
d.root = root;
d.expts = 'Functions';
d.settings = 'Settings';
d.seq = 'Sequences';
d.out = 'Outputs';
% End initialization code - DO NOT EDIT
%disp('loaded')

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


% --- Executes just before SCIn is made visible.
function SCIn_OpeningFcn(hObject, eventdata, h, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% h    structure with h and user data (see GUIDATA)
% varargin   command line arguments to SCIn (see VARARGIN)

% Choose default command line output for SCIn
h.output = hObject;
guidata(hObject, h);

% set the window position
set(groot, 'Units', 'pixels') 
screensize = get(0, 'screensize');
set(hObject, 'Units', 'pixels')
set(hObject, 'Position', [screensize(3)/2 - 200,...
                          screensize(4)/2 - 200,...
                          400, 400]);
                      
                      
disp('*** SCIn VERSION 0.1 ***');

% set global variable d: list of directories
global d
try
    rootdir
    d.root = root;
    cd(root)
catch
    d.root = pwd;
end
    
addpath(genpath(d.root))

% --- Outputs from this function are returned to the command line.
function varargout = SCIn_OutputFcn(hObject, eventdata, h) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% h    structure with h and user data (see GUIDATA)

% Get default command line output from h structure
varargout{1} = h.output;
movegui('northeast')


% --- Executes on selection change in ExptOpt.
%function ExptOpt_Callback(hObject, eventdata, h)
% hObject    handle to ExptOpt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% h    structure with h and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns ExptOpt contents as cell array
%        contents{get(hObject,'Value')} returns selected item from ExptOpt
%expts = get(hObject,'String');
%opt = get(h.ExptOpt,'Value');
%[~,h.exptFun,] = fileparts(expts{opt});
%set(h.info, 'String', 'Experiment selected');
%guidata(hObject, h)


% --- Executes during object creation, after setting all properties.
%function ExptOpt_CreateFcn(hObject, eventdata, h)
% hObject    handle to ExptOpt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% h    empty - h not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
%if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
%    set(hObject,'BackgroundColor','white');
%end
%global d
%files = dir(fullfile(d.root,d.expts,'*.m'));
%names = {files(:).name};
%set(hObject,'String',names);

% --- Executes on selection change in SettingsOpt.
function SettingsOpt_Callback(hObject, eventdata, h)
% hObject    handle to SettingsOpt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% h    structure with h and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns SettingsOpt contents as cell array
%        contents{get(hObject,'Value')} returns selected item from SettingsOpt
settings = get(hObject,'String');
opt = get(h.SettingsOpt,'Value');
[~,h.SettingsFun,~] = fileparts(settings{opt});
set(h.info, 'String', 'Settings selected');
guidata(hObject, h)

% update options
opt = 'setoptions';
eval(['h = ' h.SettingsFun '(h,opt);']);
set(h.Options,'String',h.SettingsOptions);
guidata(hObject, h)

% --- Executes during object creation, after setting all properties.
function SettingsOpt_CreateFcn(hObject, eventdata, h)
% hObject    handle to SettingsOpt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% h    empty - h not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
global d
files = dir(fullfile(d.root,d.settings,'*.m'));
names = {files(:).name};
set(hObject,'String',names);

% --- Executes on selection change in Options.
function Options_Callback(hObject, eventdata, h)
% hObject    handle to Options (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% h    structure with h and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns Options contents as cell array
%        contents{get(hObject,'Value')} returns selected item from Options
options = get(hObject,'String');
opt = get(hObject,'Value');
h.OptName = options{opt};
set(h.info, 'String', 'Option selected');
guidata(hObject, h)

% --- Executes during object creation, after setting all properties.
function Options_CreateFcn(hObject, eventdata, h)
% hObject    handle to Options (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% h    empty - h not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on button press in CreateSeq.
function CreateSeq_Callback(hObject, eventdata, h)
% hObject    handle to CreateSeq (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% h    structure with h and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of CreateSeq

% load functions if not already loaded
set(h.info, 'String', 'Loading settings...');
if ~isfield(h,'SettingsFun')
    settings = get(h.SettingsOpt,'String');
    opt = get(h.SettingsOpt,'Value');
    [~,h.SettingsFun,] = fileparts(settings{opt});
end

if ~isfield(h,'OptName')
    options = get(h.Options,'String');
    opt = get(h.Options,'Value');
    h.OptName = options{opt};
end

eval(['h = ' h.SettingsFun '(h,h.OptName);']);
set(h.info, 'String', 'Creating sequence...');
pause(0.1) % otherwise message not displayed
eval(['h = ' h.SeqFun '(h);']);
set(h.info, 'String', 'Sequence created. Save it?');
guidata(hObject, h)

% --- Executes on button press in SaveSeq.
function SaveSeq_Callback(hObject, eventdata, h)
% hObject    handle to SaveSeq (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% h    structure with h and user data (see GUIDATA)
% Hint: get(hObject,'Value') returns toggle state of SaveSeq
set(h.info, 'String', 'Saving sequence...');
filename = ['Sequence_' h.SettingsFun '_Option' h.OptName];
global d
seq = h.Seq;
settings = h.Settings;
uisave({'seq','settings'},fullfile(d.root,d.seq,filename));
set(h.info, 'String', 'Sequence saved.');

% update seqopt
files = dir(fullfile(d.root,d.seq,'*.mat'));
names = {files(:).name};
set(h.SeqOpt,'String',names);


% --- Executes on selection change in SeqOpt.
function SeqOpt_Callback(hObject, eventdata, h)
% hObject    handle to SeqOpt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% h    structure with h and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns SeqOpt contents as cell array
%        contents{get(hObject,'Value')} returns selected item from SeqOpt
seq = get(hObject,'String');
opt = get(h.SeqOpt,'Value');
h.SeqName = seq{opt};
set(h.info, 'String', 'Sequence selected');
guidata(hObject, h)

% --- Executes during object creation, after setting all properties.
function SeqOpt_CreateFcn(hObject, eventdata, h)
% hObject    handle to SeqOpt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% h    empty - h not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
global d
files = dir(fullfile(d.root,d.seq,'*.mat'));
names = {files(:).name};
set(hObject,'String',names);

% --- Executes on button press in StartStop.
function StartStop_Callback(hObject, eventdata, h)
% hObject    handle to StartStop (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% h    structure with h and user data (see GUIDATA)
global d
% get GUI handle name: necessary if 'h' is empty because not called from
% the base workspace
%GUIhname = findall(0, 'Type', 'figure', 'Tag', 'SCIn');
%h = guihandles(GUIhname);

h.d=d;

% load functions if not already loaded
set(h.info, 'String', 'Loading sequence...');
%if ~isfield(h,'exptFun')
%    expts = get(h.ExptOpt,'String');
%    opt = get(h.ExptOpt,'Value');
%    [~,h.exptFun,] = fileparts(expts{opt});
%end
h.exptFun = 'Experiment';
if ~isfield(h,'SeqName')
    seq = get(h.SeqOpt,'String');
    opt = get(h.SeqOpt,'Value');
    h.SeqName = seq{opt};
end
if ~isfield(h,'Seq') || ~isfield(h,'Settings')
    % load sequence
    A=load(fullfile(d.root,d.seq,h.SeqName));
    h.Seq = A.seq;
    h.Settings = A.settings;
end
if ~isfield(h,'startblock')
    h.startblock = '1';
end

% check if the button is pressed
if get(hObject, 'Value') == get(hObject, 'Max')

    % setup
    set(h.info, 'String', 'Setting up...');
    opt = 'setup';
    eval(['h = ' h.exptFun '(h,opt);']);
    if h.Settings.labjack
        set(h.ljhandle, 'Value', h.ljHandle);
    end
 
    
    % select blocks to run
    set(h.info, 'String', 'Setting blocks...');
    numblock = max(h.Seq.blocks);
    startblock = str2double(h.startblock);
    btrials = startblock:numblock;
    dti_all=[];
    for bt = 1:length(btrials)
        dti = find(h.Seq.blocks==btrials(bt));
        dti_all = [dti_all dti];
    end
    dti_all = sort(dti_all);
    h.Seq.signal = h.Seq.signal(:, dti_all);
    
    % save run info
    set(h.info, 'String', 'Saving run info...');
    if ~isfield(h,'subID')
        h.subID = 'test';
    end
    runinfo.seqname = h.SeqName;
    t_start = datestr(now,30);
    fname = ['runinfo_' h.subID '_startblock' h.startblock '_' t_start];
    save(fullfile(d.root,d.out,fname),'runinfo');

    % disable all GUI components except the toggle buttons
    %set(h.ExptOpt, 'Enable', 'off')
    set(h.SeqOpt, 'Enable', 'off')
    set(h.SettingsOpt, 'Enable', 'off')
    %set(h.StartStop, 'Enable', 'off')

    % startstop
    set(h.info, 'String', 'Running sequence...');
    opt = 'start';
    set(hObject, 'String', 'Stop')
    eval(['h = ' h.exptFun '(h,opt);']);
    
else

    try
        % enable all GUI components if the toggle buttons are unpressed
        if get(h.StartStop, 'Value') == get(h.StartStop, 'Min')
            %set(h.ExptOpt, 'Enable', 'on')
            set(h.SeqOpt, 'Enable', 'on')
            set(h.SettingsOpt, 'Enable', 'on')
            %set(h.StartStop, 'Enable', 'on')
        else
            %set(h.ExptOpt, 'Enable', 'off')
            set(h.SeqOpt, 'Enable', 'off')
            set(h.SettingsOpt, 'Enable', 'off')
            %set(h.StartStop, 'Enable', 'off')
        end

        % info
        set(h.info, 'String', 'Stopped.');

        % stop 
        opt = 'stop';
        eval(['h = ' h.exptFun '(h,opt);']);
    end

    % change the button text
    set(hObject, 'String', 'Start')
end
    

% --- Executes on button press in PauseResume.
function PauseResume_Callback(hObject, eventdata, h)
% hObject    handle to PauseResume (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% h    structure with h and user data (see GUIDATA)
% Hint: get(hObject,'Value') returns toggle state of PauseResume
global d
% get GUI handle name: necessary if 'h' is empty because not called from
% the base workspace
GUIhname = findall(0, 'Type', 'figure', 'Tag', 'SCIn');
h = guihandles(GUIhname);

if 0
    % load fields if they don't exist
    %if ~isfield(h,'exptFun')
    %    expts = get(h.ExptOpt,'String');
    %    opt = get(h.ExptOpt,'Value');
    %    [~,h.exptFun,] = fileparts(expts{opt});
    %end
    if ~isfield(h,'SeqName')
        seq = get(h.SeqOpt,'String');
        opt = get(h.SeqOpt,'Value');
        h.SeqName = seq{opt};
    end
    if ~isfield(h,'Seq') || ~isfield(h,'Settings')
        % load sequence
        load(fullfile(d.root,d.seq,h.SeqName));
        h.Seq = seq;
        h.Settings = settings;
    end
    if ~isfield(h,'startblock')
        h.startblock = '1';
    end
    guidata(hObject, h)
end

% check if the button is pressed
if get(hObject, 'Value') == get(hObject, 'Max')

    % info
    %try
        set(h.info, 'String', 'Paused.');
        % pause
        %opt = 'pause';
        %eval(['h = ' h.exptFun '(h,opt);']);
    %end
    

    % change the button text
    set(hObject, 'String', 'Resume')
    
else

    % info
    %try
        set(h.info, 'String', 'Running sequence...');
        % resume 
        %opt = 'resume';
        %eval(['h = ' h.exptFun '(h,opt);']);
    %end
    

    % change the button text
    set(hObject, 'String', 'Pause')
end
    
function SubID_Callback(hObject, eventdata, h)
% hObject    handle to SubID (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% h    structure with h and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of SubID as text
%        str2double(get(hObject,'String')) returns contents of SubID as a double
h.subID = get(hObject,'String');
guidata(hObject, h)


% --- Executes during object creation, after setting all properties.
function SubID_CreateFcn(hObject, eventdata, h)
% hObject    handle to SubID (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% h    empty - h not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function StartBlock_Callback(hObject, eventdata, h)
% hObject    handle to StartBlock (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% h    structure with h and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of StartBlock as text
%        str2double(get(hObject,'String')) returns contents of StartBlock as a double
h.startblock = get(hObject,'String');
guidata(hObject, h)

% --- Executes during object creation, after setting all properties.
function StartBlock_CreateFcn(hObject, eventdata, h)
% hObject    handle to StartBlock (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% h    empty - h not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on key press with focus on SCIn or any of its controls.
%function SCIn_WindowKeyPressFcn(hObject, eventdata, h)
% hObject    handle to SCIn (see GCBO)
% eventdata  structure with the following fields (see FIGURE)
%	Key: name of the key that was pressed, in lower case
%	Character: character interpretation of the key(s) that was pressed
%	Modifier: name(s) of the modifier key(s) (i.e., control, shift) pressed
% handles    structure with handles and user data (see GUIDATA)
%set(h.buttonpressed, 'String', eventdata.Key);
%set(h.buttontime, 'String', GetSecs);
%guidata(hObject, h)
%disp(['button press: ' eventdata.Key])



function vol_atten_Callback(hObject, eventdata, h)
% hObject    handle to vol_atten (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of vol_atten as text
%        str2double(get(hObject,'String')) returns contents of vol_atten as a double
h.volatten = get(hObject,'String');
guidata(hObject, h)

% --- Executes during object creation, after setting all properties.
function vol_atten_CreateFcn(hObject, eventdata, h)
% hObject    handle to vol_atten (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function entrainfreq_Callback(hObject, eventdata, h)
% hObject    handle to entrainfreq (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of entrainfreq as text
%        str2double(get(hObject,'String')) returns contents of entrainfreq as a double
h.entrainfreq = get(hObject,'String');
guidata(hObject, h)

% --- Executes during object creation, after setting all properties.
function entrainfreq_CreateFcn(hObject, eventdata, h)
% hObject    handle to entrainfreq (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
