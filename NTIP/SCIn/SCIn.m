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
dbclear if error

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
d.randseq = 'RandomLists';
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
                          screensize(4)/2 - 800,...
                          400, 600]);
                      
                      
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

if isfield(h,'Settings')
    h = rmfield(h,'Settings')
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
uisave({'seq','settings','filename'},fullfile(d.root,d.seq,filename));
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


% --- Executes on button press in CreateRand.
function CreateRand_Callback(hObject, eventdata, h)
% hObject    handle to CreateRand (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of CreateRand
global d

% CREATE
seqnames = uigetfile(fullfile(d.root,d.seq,'*.mat'),'Select sequence files','MultiSelect','on');
seqnames = seqnames(randperm(length(seqnames)));

% SAVE
set(h.info, 'String', 'Randomising sequence list...');
if ~isfield(h,'subID'); h.subID='test';end
for i = 1:length(seqnames)
    filename = ['RandSeq_' h.subID '_Seq' num2str(i) '.mat'];
    copyfile(fullfile(d.root,d.seq,seqnames{i}),fullfile(d.root,d.randseq,filename));
end
set(h.info, 'String', 'Randomisation complete.');

% update randopt
files = dir(fullfile(d.root,d.randseq,'*.mat'));
names = {files(:).name};
set(h.RandOpt,'String',names);

% --- Executes on selection change in RandOpt.
function RandOpt_Callback(hObject, eventdata, h)
% hObject    handle to RandOpt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns RandOpt contents as cell array
%        contents{get(hObject,'Value')} returns selected item from RandOpt
randseq = get(hObject,'String');
opt = get(h.RandOpt,'Value');
h.RandName = randseq{opt};
set(h.info, 'String', 'Randomised list selected');
guidata(hObject, h)

% --- Executes during object creation, after setting all properties.
function RandOpt_CreateFcn(hObject, eventdata, h)
% hObject    handle to RandOpt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
global d
files = dir(fullfile(d.root,d.randseq,'*.mat'));
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

h.exptFun = 'Experiment';
if h.mode.runseq && h.mode.runrandseq
    errordlg('first select EITHER seq OR rand-seq');
    return
elseif h.mode.runseq
    if ~isfield(h,'SeqName')
        seq = get(h.SeqOpt,'String');
        opt = get(h.SeqOpt,'Value');
    end
    A=load(fullfile(d.root,d.seq,h.SeqName));
elseif h.mode.runrandseq
    if ~isfield(h,'RandName')
        seq = get(h.RandOpt,'String');
        opt = get(h.RandOpt,'Value');
    end
    A=load(fullfile(d.root,d.randseq,h.RandName));
end

h.Seq = A.seq;
h.Settings = A.settings;
h.SeqName = A.filename;
if ~isfield(h,'startblock')
    h.startblock = '1';
end

% check if the button is pressed
if get(hObject, 'Value') == get(hObject, 'Max')

    % setup
    set(h.info, 'String', 'Setting up...');
    opt = 'setup';
    eval(['h = ' h.exptFun '(h,opt);']);
    d = h.d;
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
        prompt = 'Enter subject ID';
        h.subID = inputdlg(prompt,'Subject ID',1,{'test'});
        h.subID = h.subID{:};
    end
    %runinfo.seqname = h.SeqName;
    %t_start = datestr(now,30);
    %fname = ['runinfo_' h.subID '_startblock' h.startblock '_' t_start];
    %save(fullfile(d.root,d.out,fname),'runinfo');

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
        % stop 
        opt = 'stop';
        eval(['h = ' h.exptFun '(h,opt);']);
    end
    
    try
        % info
        set(h.info, 'String', 'Stopped.');

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
    end
    
    % change the button text
    set(hObject, 'String', 'Start')
    
    try
        global spt1
        fclose(spt1);
    end
    try
        global spt
        fclose(spt);
    end
    
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

function IntensityMean_Callback(hObject, eventdata, h)
% hObject    handle to IntensityMean (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of IntensityMean as text
%        str2double(get(hObject,'String')) returns contents of IntensityMean as a double
h.inten_mean_gui = get(hObject,'String');
guidata(hObject, h)

% --- Executes during object creation, after setting all properties.
function IntensityMean_CreateFcn(hObject, eventdata, h)
% hObject    handle to IntensityMean (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function IntensityDiff_Callback(hObject, eventdata, h)
% hObject    handle to IntensityDiff (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of IntensityDiff as text
%        str2double(get(hObject,'String')) returns contents of IntensityDiff as a double
h.inten_diff_gui = get(hObject,'String');
guidata(hObject, h)

% --- Executes during object creation, after setting all properties.
function IntensityDiff_CreateFcn(hObject, eventdata, h)
% hObject    handle to IntensityDiff (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function AudDiff_Callback(hObject, eventdata, h)
% hObject    handle to AudDiff (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of AudDiff as text
%        str2double(get(hObject,'String')) returns contents of AudDiff as a double
h.aud_diff_gui = get(hObject,'String');
guidata(hObject, h)

% --- Executes during object creation, after setting all properties.
function AudDiff_CreateFcn(hObject, eventdata, h)
% hObject    handle to AudDiff (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on button press in Setup.
function Setup_Callback(hObject, eventdata, h)
% hObject    handle to Setup (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of Setup
h.mode.setup = get(hObject,'Value');
if h.mode.setup
    set(h.SettingsOpt,'Enable','on');
    set(h.Options,'Enable','on');
    set(h.CreateSeq,'Enable','on');
    set(h.SaveSeq,'Enable','on');
    set(h.CreateRand,'Enable','on');
else
    set(h.SettingsOpt,'Enable','off');
    set(h.Options,'Enable','off');
    set(h.CreateSeq,'Enable','off');
    set(h.SaveSeq,'Enable','off');
    set(h.CreateRand,'Enable','off');
end
guidata(hObject, h)

% --- Executes on button press in runseq.
function runseq_Callback(hObject, eventdata, h)
% hObject    handle to runseq (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of runseq
h.mode.runseq = get(hObject,'Value');

if h.mode.runseq
    set(h.SeqOpt,'Enable','on');
    set(h.RandOpt,'Enable','off');
    set(h.SettingsOpt,'Enable','off');
    set(h.Options,'Enable','off');
    set(h.Setup,'value',0);
    set(h.runrandseq,'value',0);
    h.mode.runrandseq = 0;
    h.mode.setup = 0;
    set(h.CreateSeq,'Enable','off');
    set(h.SaveSeq,'Enable','off');
    set(h.CreateRand,'Enable','off');
else
    set(h.SeqOpt,'Enable','off');
end
guidata(hObject, h)

% --- Executes on button press in runrandseq.
function runrandseq_Callback(hObject, eventdata, h)
% hObject    handle to runrandseq (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of runrandseq
h.mode.runrandseq = get(hObject,'Value');
if h.mode.runrandseq
    set(h.RandOpt,'Enable','on');
    set(h.SeqOpt,'Enable','off');
    set(h.SettingsOpt,'Enable','off');
    set(h.Options,'Enable','off');
    set(h.Setup,'value',0);
    set(h.runseq,'value',0);
    h.mode.runseq = 0;
    h.mode.setup = 0;
    set(h.CreateSeq,'Enable','off');
    set(h.SaveSeq,'Enable','off');
    set(h.CreateRand,'Enable','off');
else
    set(h.RandOpt,'Enable','off');
end
guidata(hObject, h)



function DutyCycle_Callback(hObject, eventdata, h)
% hObject    handle to DutyCycle (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of DutyCycle as text
%        str2double(get(hObject,'String')) returns contents of DutyCycle as a double
h.dutycycle = get(hObject,'String');
guidata(hObject, h)

% --- Executes during object creation, after setting all properties.
function DutyCycle_CreateFcn(hObject, eventdata, h)
% hObject    handle to DutyCycle (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


