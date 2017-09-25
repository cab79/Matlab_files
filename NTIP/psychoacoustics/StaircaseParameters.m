function h = StaircaseParameters()


%H.NSUB = SUBJ NUMBER
%H.NAME = SUBJ NAME
%H.GENDER = SUBJ GENDER
%H.AGE = SUBJ AGE
%H.NOTE = OPTIONAL. SUBJ PARTICULARIT
%FILEOUT = NAME OF THE FILE WHERE THE DATA ARE SAVED
%H.NBLOCKS  =  NUMBER OF BLOCS OF THE ML PROCEDURE
%H.STIMULUSLEVEL = INIZIALIZES A VECTOR FOR THE LEVEL OF THE STIMULI
%H.SUBJCTACCURACY = INIZIALIZES A VECTOR FOR EVALUATION OF SUBJECT ACCURACY
%H.FA = INIZIALIZES A VECTOR FOR THE FALSE ALLARMS
%H.FEEDBACK= LOGICAL VALUE. DISPLAY ANSWER CORRECTNESS
%H.TEMPORARYTHRESHOLD = TEMPORARY THRESHOLDS AFTER EACH TRIAL
%H.MATSAVEDATA = INIZIALIZES A MATRIX FOR SAVING THE DATA


global h;
global d;
global myslash;
rehash toolboxcache



if ispc
    myslash = '\';
else
    myslash = '/';
end;

h.nblocks = '';
h.startinglevel = '';
h.standard= '';

h.feature='';
h.down='';
h.feedback= 0;
h.reversals = '';
h.stepsize = '';

h.fileout = '';
h.SaveResults =1;
h.tasktype = 0;
h.exppos = 1;
h.experiment = '';
h.description = '';
h.thresholdtype='';
h.feature = '';
h.reversalForthresh = '';
h.isstep = 0;
h.tasktype=0;
h.nafc=0;
h.NameStepSize='';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
mainfigure = figure (100);
set(0,'Units','pixels')
scnsize = get(0, 'ScreenSize');
hpcont = scnsize(3)/36;                     %horizontal posizion control
hsizecont = scnsize(3)/14.2;                %horizontal size control
vpcont = scnsize(4)/2;                      %vertical posizion control
vsizecont = scnsize(4)/40;                  %vertical size control
md = ([matlabroot myslash 'toolbox' myslash 'psychoacoustics' myslash 'Experiments_S']);

if exist(md, 'dir')
    d = dir([md myslash '*.m']);
else
    d = dir([cd myslash 'Experiments_S' myslash '*.m']);
end;


if ~exist('d','var')
    errordlg('Directory Experiments Staircase not found','Staircase')
    clear;
    clc;
    close;
    return;
end;

set (mainfigure, 'Position', [hpcont*7.2, vpcont*.4, (hpcont+hsizecont)*5.66 ,vpcont*1.4],...
    'Color', [0.8 0.8 0.8],...
    'NumberTitle', 'off',...
    'MenuBar', 'none',...
    'CloseRequestFcn',{@closeGUImasolo},...
    'Name', ' Staircase');

f = uimenu('Label','&Select Experiment');

for i = 1:size(d,1)
    fun = ['{@exp' ',' num2str(i)  '}'];
    uimenu(f,'Label',d(i).name,'Callback',eval(fun));
end;

f2 = uimenu('Label','&Edit Experiment');

for i = 1:size(d,1)
    fun = ['{@editexp' ',' num2str(i)  '}'];
    uimenu(f2,'Label',d(i).name,'Callback',eval(fun));
end;

uicontrol(mainfigure, 'Style', 'text',...
    'String', 'Experiment',...
    'Position',[hpcont vpcont+140  hpcont+hsizecont*6 vpcont-(vpcont-vsizecont)]);

uicontrol(mainfigure, 'Style', 'text',...
    'String', h.experiment,...
    'FontSize',10,'FontWeight','bold',...
    'Position',[hpcont vpcont+131  hpcont+hsizecont*6 15],...
    'Tag', 'ExpLabel');

uicontrol(mainfigure, 'Style', 'edit',...
    'String', h.description,...
    'Max',2,'Min',0,...
    'HorizontalAlignment','left',...
    'Position',[hpcont vpcont+50  hpcont+hsizecont*7 80],...
    'Tag', 'description');

uicontrol(mainfigure, 'Style', 'text',...
    'Position',[hpcont vpcont  hpcont+hsizecont vpcont-(vpcont-vsizecont)],...
    'String', 'number');

uicontrol(mainfigure, 'Style', 'edit',...
    'String','',...
    'Position',[hpcont vpcont-vsizecont  hpcont+hsizecont vpcont-(vpcont-vsizecont)],...
    'BackgroundColor', [1 1 1],...
    'Tag', 'nsub');

uicontrol(mainfigure, 'Style','text',...
    'Position',[hpcont+hsizecont*2 vpcont  hpcont+hsizecont vpcont-(vpcont-vsizecont)],...
    'String', 'name');

uicontrol(mainfigure, 'Style', 'edit',...
    'String', '',...
    'Position',[hpcont+hsizecont*2 vpcont-vsizecont  hpcont+hsizecont vpcont-(vpcont-vsizecont)],...
    'BackgroundColor', [1 1 1],...
    'Tag', 'name');

uicontrol(mainfigure, 'Style','text',...
    'Position',[hpcont+hsizecont*4 vpcont  hpcont+hsizecont vpcont-(vpcont-vsizecont)],...
    'String', 'sex');

uicontrol(mainfigure, 'Style', 'edit',...
    'String', '',...
    'Position',[hpcont+hsizecont*4 vpcont-vsizecont  hpcont+hsizecont vpcont-(vpcont-vsizecont)],...
    'BackgroundColor', [1 1 1],...
    'Tag', 'gender');

uicontrol(mainfigure, 'Style','text',...
    'Position',[hpcont+hsizecont*6 vpcont  hpcont+hsizecont vpcont-(vpcont-vsizecont)],...
    'String', 'age');

uicontrol(mainfigure, 'Style', 'edit',...
    'String', '',...
    'Position',[hpcont+hsizecont*6 vpcont-vsizecont  hpcont+hsizecont vpcont-(vpcont-vsizecont)],...
    'BackgroundColor', [1 1 1],...
    'Tag', 'age');

uicontrol(mainfigure, 'Style','text',...
    'Position',[hpcont vpcont-vsizecont*3.5  hpcont+hsizecont/10 vpcont-(vpcont-vsizecont)],...
    'String', 'note');

uicontrol(mainfigure, 'Style', 'edit',...
    'String', '',...
    'Position',[hpcont+hsizecont/2 vpcont-vsizecont*4  hpcont+hsizecont*2.5 vpcont-(vpcont-vsizecont*2)],...
    'BackgroundColor', [1 1 1],...
    'Tag', 'note');

uicontrol(mainfigure, 'Style','text',...
    'Position',[hpcont+hsizecont*3.5 vpcont-vsizecont*3.5  hpcont+hsizecont vpcont-(vpcont-vsizecont)],...
    'String', 'File data');

uicontrol(mainfigure, 'Style', 'edit',...
    'String', h.fileout,...
    'HorizontalAlignment','Center',...
    'Position',[hpcont+hsizecont*4.5 vpcont-vsizecont*3.5  hpcont+hsizecont vpcont-(vpcont-vsizecont)],...
    'BackgroundColor', [1 1 1],...
    'Tag', 'fileout');

uicontrol('Style','checkbox','String','Save Results',...
    'pos',[hpcont+hsizecont*6.3 vpcont-vsizecont*3.5  hpcont+hsizecont/1.5 vpcont-(vpcont-vsizecont)],...
    'Tag', 'SaveResults' ,'Value',1,'CallBack', {@SaveResults_Callback});

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



uicontrol(mainfigure, 'Style','text',...
    'Position',[hpcont+hsizecont vpcont-vsizecont*7  hpcont+hsizecont/3 vpcont-(vpcont-vsizecont)],...
    'String', 'n. of blocks');

uicontrol(mainfigure, 'Style', 'edit',...
    'String', h.nblocks,...
    'Position',[hpcont+hsizecont vpcont-vsizecont*8  hpcont+hsizecont/3 vpcont-(vpcont-vsizecont)],...
    'BackgroundColor', [1 1 1],...
    'Tag', 'nblocks');

uicontrol(mainfigure, 'Style', 'text',...
    'String', 'starting level',...
    'Position',[hpcont+hsizecont*2 vpcont-vsizecont*7  hpcont+hsizecont/3 vpcont-(vpcont-vsizecont)]);

uicontrol(mainfigure, 'Style', 'edit',...
    'String', h.startinglevel,...
    'BackgroundColor', [1 1 1],...
    'Position',[hpcont+hsizecont*2 vpcont-vsizecont*8  hpcont+hsizecont/3 vpcont-(vpcont-vsizecont)],...
    'Tag', 'startinglevel');

uicontrol(mainfigure, 'Style', 'text',...
    'String', 'standard level',...
    'Position',[hpcont+hsizecont*3 vpcont-vsizecont*7  hpcont+hsizecont/2 vpcont-(vpcont-vsizecont)]);

uicontrol(mainfigure, 'Style', 'edit',...
    'String', h.standard,...
    'BackgroundColor', [1 1 1],...
    'Position',[hpcont+hsizecont*3 vpcont-vsizecont*8  hpcont+hsizecont/3 vpcont-(vpcont-vsizecont)],...
    'Tag', 'standard');

buttongroup = uibuttongroup('visible','off');

yes_nobutton = uicontrol('Style','Radio','String','yes/no',...
    'pos',[hpcont+hsizecont*4.3 vpcont-vsizecont*7  hpcont+hsizecont/1.5 vpcont-(vpcont-vsizecont)],...
    'Tag', 'yes_nobutton',...
    'parent',buttongroup,'HandleVisibility','off');

nafcbutton = uicontrol('Style','Radio','String','nAFC',...
    'pos',[hpcont+hsizecont*4.3 vpcont-vsizecont*8  hpcont+hsizecont/1.5 vpcont-(vpcont-vsizecont)],...
    'Tag', 'nafcbutton',...
    'parent',buttongroup,'HandleVisibility','off');

set(buttongroup,'SelectionChangeFcn',@selcbk);
set(buttongroup,'Visible','on');




uicontrol(mainfigure, 'Style','text',...
    'Position',[hpcont+hsizecont*5 vpcont-vsizecont*7  hpcont+hsizecont/1.5 vpcont-(vpcont-vsizecont)],...
    'String', 'nAFC','Tag','textnafc','visible','off');

uicontrol(mainfigure, 'Style', 'edit',...
    'String', h.nafc,...
    'Position',[hpcont+hsizecont*5 vpcont-vsizecont*8  hpcont+hsizecont/1.5 vpcont-(vpcont-vsizecont)],...
    'BackgroundColor', [1 1 1],...
    'Tag', 'nafc','visible','off');



uicontrol('Style','checkbox','String','Feedback',...
    'pos',[hpcont+hsizecont*6.2 vpcont-vsizecont*6.4  hpcont+hsizecont/1.5 vpcont-(vpcont-vsizecont)],...
    'Tag', 'feedback' ,'Value',h.feedback,'HandleVisibility','on');

uicontrol(mainfigure, 'Style','text',...
    'Position',[hpcont+hsizecont*1.9 vpcont-vsizecont*11 hpcont+hsizecont*2 vpcont-(vpcont-vsizecont)],...
    'String', 'Reversals','Tag','Namereversals','visible','on');

uicontrol(mainfigure, 'Style','edit',...
    'BackgroundColor', [1 1 1],...
    'String', h.reversals,...
    'Position',[hpcont+hsizecont*1.9 vpcont-vsizecont*11.7 hpcont+hsizecont*2 vpcont-(vpcont-vsizecont)],...
    'Tag','reversals','HorizontalAlignment','left','visible','on');


buttongroup = uibuttongroup('visible','off');

step = uicontrol('Style','Radio','String','step',...
    'pos',[hpcont+hsizecont*1.9 vpcont-vsizecont*13.8  hpcont+hsizecont vpcont-(vpcont-vsizecont)],...
    'Tag', 'stepbutton',...
    'parent',buttongroup,'HandleVisibility','off');

factor = uicontrol('Style','Radio','String','factor',...
    'pos',[hpcont+hsizecont*1.9 vpcont-vsizecont*14.8  hpcont+hsizecont vpcont-(vpcont-vsizecont)],...
    'Tag', 'factorbutton',...
    'parent',buttongroup,'HandleVisibility','off');


set(buttongroup,'SelectionChangeFcn',@stepfactor);
set(buttongroup,'Visible','on');

uicontrol(mainfigure, 'Style','text',...
    'Position',[hpcont+hsizecont*2.6 vpcont-vsizecont*13.8 hpcont+hsizecont*1.3 vpcont-(vpcont-vsizecont)],...
    'String', h.NameStepSize,'Tag','NameStepsize','visible','on');


uicontrol(mainfigure, 'Style','edit',...
    'BackgroundColor', [1 1 1],...
    'String',h.stepsize,...
    'Position',[hpcont+hsizecont*2.6 vpcont-vsizecont*14.5 hpcont+hsizecont*1.3 vpcont-(vpcont-vsizecont)],...
    'Tag','stepsize','HorizontalAlignment','left','visible','on');


buttongroup = uibuttongroup('visible','off');

Arithmetic = uicontrol('Style','Radio','String','Arithmetic',...
    'pos',[hpcont+hsizecont*4.8 vpcont-vsizecont*11  hpcont+hsizecont/1.5 vpcont-(vpcont-vsizecont)],...
    'Tag', 'Arithmetic',...
    'parent',buttongroup,'HandleVisibility','off');

Geometric = uicontrol('Style','Radio','String','Geometric',...
    'pos',[hpcont+hsizecont*4.8 vpcont-vsizecont*12  hpcont+hsizecont/1.4 vpcont-(vpcont-vsizecont)],...
    'Tag', 'Geometric',...
    'parent',buttongroup);
Median = uicontrol('Style','Radio','String','Median',...
    'pos',[hpcont+hsizecont*4.8 vpcont-vsizecont*13  hpcont+hsizecont/1.5 vpcont-(vpcont-vsizecont)],...
    'Tag', 'Median',...
    'parent',buttongroup);



set(buttongroup,'SelectionChangeFcn',@selthreshold);
set(buttongroup,'Visible','on');

buttongroup = uibuttongroup('visible','off');
updown = uicontrol('Style','Radio','String','SimpleUpdown',...
    'pos',[hpcont vpcont-vsizecont*12  hpcont+hsizecont/1.5 vpcont-(vpcont-vsizecont)],...
    'Tag', 'SimpleUpdown',...
    'parent',buttongroup,'HandleVisibility','off');


limits = uicontrol('Style','Radio','String','MethodOfLimits',...
    'pos',[hpcont vpcont-vsizecont*11  hpcont+hsizecont/1.2 vpcont-(vpcont-vsizecont)],...
    'Tag', 'MethodOfLimits',...
    'parent',buttongroup,'HandleVisibility','off');

twodown = uicontrol('Style','Radio','String','TwoDownOneUp',...
    'pos',[hpcont vpcont-vsizecont*13  hpcont+hsizecont/1.2 vpcont-(vpcont-vsizecont)],...
    'Tag', 'TwoDownOneUp',...
    'parent',buttongroup,'HandleVisibility','off');

threedown = uicontrol('Style','Radio','String','ThreeDownOneUp',...
    'pos',[hpcont vpcont-vsizecont*14  hpcont+hsizecont/1.2 vpcont-(vpcont-vsizecont)],...
    'Tag', 'ThreeDownOneUp',...
    'parent',buttongroup,'HandleVisibility','off');
fourdown = uicontrol('Style','Radio','String','FourDownOneUp',...
    'pos',[hpcont vpcont-vsizecont*15  hpcont+hsizecont/1.2 vpcont-(vpcont-vsizecont)],...
    'Tag', 'FourDownOneUp',...
    'parent',buttongroup,'HandleVisibility','off');

set(buttongroup,'SelectionChangeFcn',@selmethod);
set(buttongroup,'Visible','on');


uicontrol(mainfigure, 'Style','text',...
    'Position',[hpcont+hsizecont*6 vpcont-vsizecont*11 hpcont+hsizecont vpcont-(vpcont-vsizecont)],...
    'String', 'Reversal for threshold','Tag','Namereversals','visible','on');

uicontrol(mainfigure, 'Style','edit',...
    'BackgroundColor', [1 1 1],...
    'String', h.reversalForthresh,...
    'Position',[hpcont+hsizecont*6 vpcont-vsizecont*12 hpcont+hsizecont/1.5 vpcont-(vpcont-vsizecont)],...
    'Tag','reversalForthresh','visible','on');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
uipanel('Title','Demographic data',...
    'BackgroundColor',[.9 .9 .9],...
    'Position', [.01 .56 .98 .22],...
    'FontSize',8,'FontAngle','italic');

uipanel('Title','Experiment features','FontSize',8,...
    'BackgroundColor',[.9 .9 .9],...
    'Position', [.01 .4 .98 .15],...
    'FontSize',8,'FontAngle','italic');

uipanel('Title','Staircase features','FontSize',8,...
    'BackgroundColor',[.9 .9 .9],...
    'Position', [.01 .15 .98 .25],...
    'FontSize',8,'FontAngle','italic');

uicontrol(mainfigure, 'Style', 'pushbutton',...
    'Position',[hpcont vpcont-vsizecont*18  hpcont+hsizecont vpcont-(vpcont-vsizecont*2)],...
    'String', 'START',...
    'FontSize',10,'FontWeight','bold',...
    'BackgroundColor', [.9 .9 .9],...
    'CallBack', {@START_Callback});

uicontrol(mainfigure, 'Style', 'pushbutton',...
    'Position',[hpcont+hsizecont*4 vpcont-vsizecont*18  hpcont+hsizecont+40 vpcont-(vpcont-vsizecont*2)],...
    'String', 'SAVE DEFAULTS',...
    'FontSize',10,'FontWeight','bold',...
    'BackgroundColor', [.9 .9 .9],...
    'CallBack', {@ChangeDefault_Callback});

uicontrol(mainfigure, 'Style', 'pushbutton',...
    'Position',[hpcont+hsizecont*6 vpcont-vsizecont*18  hpcont+hsizecont vpcont-(vpcont-vsizecont*2)],...
    'String', 'Cancel',...
    'FontSize',10,'FontAngle','italic',...
    'BackgroundColor', [.9 .9 .9],...
    'CallBack', {@Cancel_Callback});
exp(0,0);
uiwait;


function stepfactor(source,eventdata)
global h
mh = guihandles(gcf);
sel= get(eventdata.NewValue,'String');
switch sel
    case 'step'
        set(mh.NameStepsize,'string','Step Size');
        h.isstep=1;
    case 'factor'
        set(mh.NameStepsize,'string','Factor');
        h.isstep=0;
end

function selcbk(source,eventdata)
mh = guihandles(gcf);
sel= get(eventdata.NewValue,'String');
switch sel
    case 'yes/no'
        set(mh.nafc,'Visible', 'off');
        set(mh.textnafc,'Visible','off');
        set(mh.nafc,'String', 0);
        set(mh.feedback,'Value', 0);

    case 'nAFC'
        set(mh.nafc,'Visible', 'on');
        set(mh.textnafc,'Visible','on');
        set(mh.nafc,'String', 2);
        set(mh.feedback,'Value', 1);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function selmethod(source,eventdata)

mh = guihandles(gcf);
sel= get(eventdata.NewValue,'String');

switch sel
    case 'MethodOfLimits'
        down=1;
        set(mh.reversals,'Style', 'text');
        set(mh.reversals,'BackgroundColor', [0.8 0.8 0.8]);
        set(mh.reversals,'String', '1');
        set(mh.reversalForthresh,'Style', 'text');
        set(mh.reversalForthresh,'BackgroundColor', [0.8 0.8 0.8]);
        set(mh.reversalForthresh,'String', '1');
        set(mh.Arithmetic,'value',1');
        set(mh.Median,'visible','off');
        set(mh.Geometric,'visible','off');
        set(mh.stepsize,'String', ' ');

    case 'SimpleUpdown'
        down=1;
        set(mh.reversals,'Style', 'edit');
        set(mh.reversals,'BackgroundColor', [1 1 1]);
        set(mh.reversalForthresh,'Style', 'edit');
        set(mh.reversalForthresh,'BackgroundColor', [1 1 1]);
        set(mh.Median,'visible','on');
        set(mh.Geometric,'visible','on');

    case 'TwoDownOneUp'
        down=2;
        set(mh.reversals,'Style', 'edit');
        set(mh.reversals,'BackgroundColor', [1 1 1]);
        set(mh.reversalForthresh,'Style', 'edit');
        set(mh.reversalForthresh,'BackgroundColor', [1 1 1]);
        set(mh.Median,'visible','on');
        set(mh.Geometric,'visible','on');

    case 'ThreeDownOneUp'
        down=3;
        set(mh.reversals,'Style', 'edit');
        set(mh.reversals,'BackgroundColor', [1 1 1]);
        set(mh.reversalForthresh,'Style', 'edit');
        set(mh.reversalForthresh,'BackgroundColor', [1 1 1]);
        set(mh.Median,'visible','on');
        set(mh.Geometric,'visible','on');

    case 'FourDownOneUp'
        down=4;
        set(mh.reversals,'Style', 'edit');
        set(mh.reversals,'BackgroundColor', [1 1 1]);
        set(mh.reversalForthresh,'Style', 'edit');
        set(mh.reversalForthresh,'BackgroundColor', [1 1 1]);
        set(mh.Median,'visible','on');
        set(mh.Geometric,'visible','on');
end
global h;
h.feature=sel;
h.down=down;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function selthreshold(source,eventdata)
global h;
mh = guihandles(gcf);
sel= get(eventdata.NewValue,'String');
h.thresholdtype=sel;

%%%%%%%%%%%%%%%%%%
function ChangeDefault_Callback(a,b)
global h;
savestart (1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function START_Callback(a,b)
global h;
savestart(0);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Cancel_Callback(a,b)
global h;
if gcf
    delete (gcf);
end;
h=[];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function savestart(savedata)
global h;
global myslash;
mh = guihandles(gcf);

h.nsub = str2double(get(mh.nsub,'String'));
h.name = get(mh.name,'String');
h.gender = get(mh.gender,'String');
h.age = str2double( get(mh.age,'String'));
h.note = get(mh.note,'String');
h.fileout = get(mh.fileout,'String');
h.nblocks = str2double(get(mh.nblocks,'String'));
h.startinglevel = str2double(get(mh.startinglevel,'String'));
h.standard = str2double(get(mh.standard,'String'));
h.reversals = get(mh.reversals,'String');
h.stepsize = get(mh.stepsize,'String');
% h.isstep = str2double(get(mh.isstep,'String'));
h.reversals = returnnum(h.reversals);
h.stepsize = returnnum(h.stepsize);
h.reversalForthresh = str2double(get(mh.reversalForthresh,'String'));

if get(mh.FourDownOneUp,'Value')
    h.feature= 'FourDownOneUp';
    h.down=4;
end

if get(mh.ThreeDownOneUp,'Value')
    h.feature= 'ThreeDownOneUp';
    h.down=3;
end
if get(mh.TwoDownOneUp,'Value')
    h.feature= 'TwoDownOneUp';
    h.down=2;
end

if get(mh.SimpleUpdown,'Value')
    h.feature= 'SimpleUpdown';
    h.down=1;
end


if get(mh.MethodOfLimits,'Value')
    h.feature= 'MethodOfLimits';
    h.down=1;
end




if get(mh.Arithmetic,'Value')
    h.thresholdtype ='Arithmetic';
end

if get(mh.Geometric,'Value')
    h.thresholdtype ='Geometric';
end

if get(mh.Median,'Value')
    h.thresholdtype ='Median';
end

h.nafc =str2double( get(mh.nafc,'String'));
disp (h.nafc)
if h.nafc
    h.tasktype=1;
else
    h.tasktype=0;
end

if iscell(get(mh.feedback,'Value'))
    tempVal=get(mh.feedback,'Value');
    tempMax=get(mh.feedback,'Max');
    h.feedback= tempVal{1}==tempMax{1};

else
    h.feedback = get(mh.feedback,'Value')==get(mh.feedback,'Max');
end

h.description = get(mh.description,'String');


if length(h.reversals) ~= length(h.stepsize)
    errordlg('The number of reversals and the number of step size must be the same','Staircase')
    return
end

if h.reversalForthresh > sum(h.reversals)
    errordlg('Reversal for threshold cannot be bigger than the sum of reversals','Staircase')
    return
end

h.NameStepSize =get(mh.NameStepsize,'string');
if strcmp( h.NameStepSize ,'Step Size')
    h.isstep=1;
else
    h.isstep=0;
end



if savedata
    mddefaults = ([matlabroot myslash 'toolbox' myslash 'psychoacoustics' myslash 'Defaults_S']);
    if exist(mddefaults, 'dir')
        defaults = [mddefaults ,myslash];
    else
        defaults = [cd ,myslash 'Defaults_S' myslash];
    end
    filexppar = [defaults , h.experiment(1:end-2),'Par_S.m'];
    if  exist(filexppar,'file')
        fileexist=1;
        answ = questdlg('Default parameters for this experiment already exists. Do you want to replace it?', 'Staircase','Yes','No','Yes');
    else
        fileexist=0;
    end
    if ~fileexist || strcmp(answ, 'Yes')
        fid3 =    fopen(filexppar, 'w');
        filepar = [h.experiment(1:end-2), 'Par_S(mh)'];
        fprintf(fid3, '%s\t','function ');
        fprintf(fid3, '%s\n',filepar);        fprintf(fid3, '%s\t','set(mh.nblocks,''String'',');
        fprintf(fid3, '%g\t',h.nblocks);
        fprintf(fid3, '%s\n',');');
        fprintf(fid3, '%s','set(mh.fileout,''String'',''');
        fprintf(fid3, '%s',h.fileout);
        fprintf(fid3, '%s\n',''');');
        fprintf(fid3, '%s\t','set(mh.startinglevel,''String'',');
        fprintf(fid3, '%g\t',h.startinglevel);
        fprintf(fid3, '%s\n',');') ;
        fprintf(fid3, '%s\t','set(mh.standard,''String'',');
        fprintf(fid3, '%g\t',h.standard);
        fprintf(fid3, '%s\n',');') ;
        fprintf(fid3, '%s\t','set(mh.reversals,''String'',');
        fprintf(fid3, '%s','''') ;
        fprintf(fid3, '%g ',h.reversals);
        fprintf(fid3, '%s','''') ;
        fprintf(fid3, '%s\n',');') ;
        fprintf(fid3, '%s\t','set(mh.stepsize,''String'',');
        fprintf(fid3, '%s','''') ;
        fprintf(fid3, '%g ',h.stepsize);
        fprintf(fid3, '%s','''') ;
        fprintf(fid3, '%s\n',');') ;
        fprintf(fid3, 'set(mh.NameStepsize,''String'', ''');
        fprintf(fid3, '%s',h.NameStepSize);
        fprintf(fid3, ''');\n') ;




        if h.isstep

            fprintf(fid3, '%s\t','set(mh.stepbutton,''value'',1');
        else
            fprintf(fid3, '%s\t','set(mh.factorbutton,''value'',1');
        end
  fprintf(fid3, '%s\n',');');



        switch h.feature
            case 'MethodOfLimits'
                fprintf(fid3, '%s\t','set(mh.MethodOfLimits,''value'',1');
                fprintf(fid3, '%s\n',');');


                fprintf(fid3, '%s\t',   'set(mh.reversals,''Style'',''text''');
                fprintf(fid3, '%s\n',');');
                fprintf(fid3, '%s\t',   'set(mh.reversals,''BackgroundColor'', [0.8 0.8 0.8]');
                fprintf(fid3, '%s\n',');');
                fprintf(fid3, '%s\t',   'set(mh.reversals,''String'', 1');
                fprintf(fid3, '%s\n',');');
                fprintf(fid3, '%s\t',   'set(mh.reversalForthresh,''Style'',''text''');
                fprintf(fid3, '%s\n',');');
                fprintf(fid3, '%s\t',   'set(mh.reversalForthresh,''BackgroundColor'',   [0.8 0.8 0.8]');
                fprintf(fid3, '%s\n',');');
                fprintf(fid3, '%s\t',   'set(mh.reversalForthresh,''String'', 1');
                fprintf(fid3, '%s\n',');');
                fprintf(fid3, '%s\t','set(mh.Arithmetic,''value'',1');
                fprintf(fid3, '%s\n',');');
                fprintf(fid3, '%s\t','set(mh.Median,''visible'',''off''');
                fprintf(fid3, '%s\n',');');
                fprintf(fid3, '%s\t','set(mh.Geometric,''visible'',''off''');
                fprintf(fid3, '%s\n',');');
            case 'SimpleUpdown'
                fprintf(fid3, '%s\t','set(mh.SimpleUpdown,''value'',1');
                fprintf(fid3, '%s\n',');');
            case 'TwoDownOneUp'
                fprintf(fid3, '%s\t','set(mh.TwoDownOneUp,''value'',1');
                fprintf(fid3, '%s\n',');');
            case 'ThreeDownOneUp'
                fprintf(fid3, '%s\t','set(mh.ThreeDownOneUp,''value'',1');
                fprintf(fid3, '%s\n',');');
            case 'FourDownOneUp'
                fprintf(fid3, '%s\t','set(mh.FourDownOneUp,''value'',1');
                fprintf(fid3, '%s\n',');');
        end

        switch h.thresholdtype
            case 'Arithmetic'
                fprintf(fid3, '%s\t','set(mh.Arithmetic,''value'',1');
                fprintf(fid3, '%s\n',');');
            case 'Geometric'
                fprintf(fid3, '%s\t','set(mh.Geometric,''value'',1');
                fprintf(fid3, '%s\n',');');
            case 'Median'
                fprintf(fid3, '%s\t','set(mh.Median,''value'',1');
                fprintf(fid3, '%s\n',');');
        end

        if h.nafc
            fprintf(fid3, '%s\n','set(mh.nafcbutton, ''value'' ,1); ');
            fprintf(fid3, '%s\n','set(mh.textnafc,''Visible'',''on''); ');
            fprintf(fid3, '%s\n','set(mh.nafc,''Visible'', ''on''); ');

            fprintf(fid3, '%s\t','set(mh.nafc, ''String'' , ');
            fprintf(fid3, '%g\t',h.nafc);
            fprintf(fid3, '%s\n',');')  ;
        else
            fprintf(fid3, '%s\n','set(mh.yes_nobutton, ''value'' ,1); ');
            fprintf(fid3, '%s\n','set(mh.nafc, ''String'' , 0);')  ;
        end

        fprintf(fid3, '%s\t','set(mh.reversalForthresh,''String'',');
        fprintf(fid3, '%g\t',h.reversalForthresh);
        fprintf(fid3, '%s\n',');') ;
        fprintf(fid3, '%s\t','set(mh.feedback, ''value'' , ');
        fprintf(fid3, '%g\t',h.feedback);
        fprintf(fid3, '%s\n',');');
        fprintf(fid3, '%s\t','set(mh.SaveResults, ''value'' , ');
        fprintf(fid3, '%g\t',h.SaveResults);
        fprintf(fid3, '%s\n',');');
        fclose(fid3);
    end

    return
else

end;

delete (gcf);
if gcf
    delete (gcf);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function closeGUImasolo(src,evnt)
global h;
selection = questdlg('Do you want to close Staircase?',...
    'Staircase',...
    'Yes','No','Yes');
switch selection,
    case 'Yes'
        h=[];
        delete(gcf)
    case 'No'
        return
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function SaveResults_Callback(hObject, eventdata, handles)
global h;
if (get(hObject,'Value') == get(hObject,'Max'))
    h.SaveResults = 1;
else
    h.SaveResults = 0;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function exp(source,eventdata,indx)
global h;
global d;
global myslash;
mddefaults = ([matlabroot myslash 'toolbox' myslash 'psychoacoustics' myslash 'Defaults_S']);


if exist(mddefaults, 'dir')
    defaults = [mddefaults ,myslash];
else
    defaults = [cd ,myslash 'Defaults_S' myslash];
end
filexp = [ defaults 'lastexp.txt'];


if nargin<3

    fid =    fopen(filexp, 'r');
    if fid>-1
        indx = fscanf(fid, '%g');
    else
        indx =1;
    end

    if exist(filexp, 'file')
        fclose(fid);
    else
        warndlg('You are running the procedure outside the psychoacustics directory. Either enter the parameters or change directory and try again','Staircase procedure');
    end
end

mh = guihandles(gcf);
if size(d,1)< indx
    indx=1;
    msgbox('Last experiment defaults have probably been deleted. The defaults refer to the first experiment in the list.','Staircase','Warn')
end;

filesel = d(indx).name ;
h.experiment = filesel(1:length(filesel)-2) ;
set(mh.ExpLabel,'String', h.experiment);
descr = help(h.experiment);
set(mh.description,'String', descr);
funpar = [h.experiment(1:end-2) 'Par_S(mh)'];
filepar= [h.experiment(1:end-2) 'Par_S.m'];

if  exist(filepar, 'file')

    eval(funpar);

else
    set(mh.nblocks, 'String', '');
    set(mh.standard, 'String', '');
    set(mh.feedback   , 'Value', 0);
    set(mh.fileout   , 'String', '');
    set(mh.SaveResults   ,'Value', 0);
    set(mh.MethodOfLimits,'value', 1);
    set(mh.Median,'value', 1);
    set(mh.reversals,'String', '');
    set(mh.stepsize,'String', '');
    warndlg(['File ' filepar ' doesn''t exist. Press SAVE DAFAULTS if you want to save defaults for this experiment'], 'Staircase');
end

fid2 =    fopen(filexp, 'w');
fprintf(fid2, '%g', indx);
fclose(fid2);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function editexp(source,eventdata,indx)
global d;

filesel = d(indx).name;
edit (filesel);

