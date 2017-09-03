function h = MLParameters()

%MLPARAMETERS Computes the parameters for Green's ML procedure.


%H.NSUB = SUBJ NUMBER
%H.NAME = SUBJ NAME
%H.GENDER = SUBJ GENDER
%H.AGE = SUBJ AGE
%H.NOTE = OPTIONAL. SUBJ PARTICULARIT
%FILEOUT = NAME OF THE FILE WHERE THE DATA ARE SAVED
%H.NTRIALS = NUMBER OF TRIALS OF THE ML PROCEDURE
%H.NBLOCKS  =  NUMBER OF BLOCS OF THE ML PROCEDURE
%H.NAFC =  NUMBER OF ALTERNATIVE FOR CHOICE
%H.STARTINGLEVEL = INTENSITY OF THE FIRST STIMULUS
%H.BETA = PSYCHOMETRIC FUNCTION SLOPE
%H.FIRSTMIDPOINT = PSYCHOMETRIC FUNCTION FIRST MIDPOINT
%H.LASTMIDPOINT = PSYCHOMETRIC FUNCTION  LAST MIDPOINT
%H.HYPNR = NUMBER OF HYPOTHESES GENERATED FOR THE THRESHOLD ESTIMATION
%H.LAMBDA = EXPECTED PROPORTION OF ATTENTIONAL LAPSES
%H.LOGSCALE = LOGICAL VALUE. HYPOTHESES LOGARITHMICALLY SPACED =1 OR NOT =0
%H.SWEETPOINT = LOGICAL VALUE. SWEETPOINT NEEDED =1 OR NOT =0
%H.P_TARGET = TARGET PROBABILITY
%H.GAMMA = THE LOWER LIMIT OF THE PSYCHOMETRIC FUNCTION
%H.CATCHTRIAL = PROPORTION OF CATCH TRIAL TO PRESENT DURING THE EXPERIMENT
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


h.ntrials = '';
h.nblocks = '';
h.startinglevel = '';
h.standard= '';
h.catchtrial = '';
h.nafc = 0;
h.p_target = '';
h.beta = '';
h.firstmidpoint = '';
h.lastmidpoint = '';
h.hypnr = '';
h.islog = 0;
h.feedback= 0;
h.repeatft= 0;
h.lambda = 0;
h.fileout = '';
h.SaveResults =1;
h.tasktype = 0;
h.swp= 0;
h.functiontype='logistic';
h.logscale = 0;
h.exppos = 1;
h.experiment = '';
h.description = '';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
mainfigure = figure (100);
set(0,'Units','pixels')
scnsize = get(0, 'ScreenSize');
hpcont = scnsize(3)/36;                     %horizontal posizion control
hsizecont = scnsize(3)/14.2;                %horizontal size control
vpcont = scnsize(4)/2;                      %vertical posizion control
vsizecont = scnsize(4)/40;                  %vertical size control
md = ([matlabroot myslash 'toolbox' myslash 'psychoacoustics' myslash 'Experiments_MLP']);

if exist(md, 'dir')
    d = dir([md myslash '*.m']);
else
    d = dir([cd myslash 'Experiments_MLP'  myslash '*.m']);
end;
if ~exist('d','var')
    errordlg('Directory Experiments not found','MLP Procedure')
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
    'Name', ' Maximum Likelihood procedure');

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

uicontrol(mainfigure, 'Style','text',...
    'Position',[hpcont vpcont-vsizecont*7  hpcont+hsizecont/2 vpcont-(vpcont-vsizecont)],...
    'String', 'trials per block');

uicontrol(mainfigure, 'Style', 'edit',...
    'String', h.ntrials,...
    'Position',[hpcont vpcont-vsizecont*8  hpcont+hsizecont/3 vpcont-(vpcont-vsizecont)],...
    'BackgroundColor', [1 1 1],...
    'Tag', 'ntrials');

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
    'String', 'standard (or catch) level)',...
    'Position',[hpcont+hsizecont*3 vpcont-vsizecont*7  hpcont+hsizecont/2.8 vpcont-(vpcont-vsizecont)]);

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

txtcatch = uicontrol(mainfigure, 'Style','text',...
    'Position',[hpcont+hsizecont*5 vpcont-vsizecont*7  hpcont+hsizecont vpcont-(vpcont-vsizecont)],...
    'String', 'catch prop','Tag','textcatchtrial','visible','off');

uicontrol(mainfigure, 'Style', 'edit',...
    'String', h.catchtrial,...
    'Position',[hpcont+hsizecont*5.1 vpcont-vsizecont*8  hpcont+hsizecont/1.5 vpcont-(vpcont-vsizecont)],...
    'BackgroundColor', [1 1 1],...
    'Tag', 'catchprop','visible','off');

uicontrol(mainfigure, 'Style','text',...
    'Position',[hpcont+hsizecont*5 vpcont-vsizecont*7  hpcont+hsizecont/1.5 vpcont-(vpcont-vsizecont)],...
    'String', 'nAFC','Tag','textnafc','visible','off');

uicontrol(mainfigure, 'Style', 'edit',...
    'String', h.nafc,...
    'Position',[hpcont+hsizecont*5 vpcont-vsizecont*8  hpcont+hsizecont/1.5 vpcont-(vpcont-vsizecont)],...
    'BackgroundColor', [1 1 1],...
    'Tag', 'nafc','visible','off');

uicontrol(mainfigure, 'Style', 'pushbutton',...
    'Position',[hpcont+hsizecont*6.2 vpcont-vsizecont*8  hpcont+hsizecont/1.5 vpcont-(vpcont-vsizecont)],...
    'String', 'Update ptarget',...
    'FontSize',10,...
    'BackgroundColor', [.9 .9 .9],...
    'visible','off',...
    'Tag', 'ptargetupdate',...
    'CallBack', {@Updateptarget_Callback});

uicontrol('Style','checkbox','String','Feedback',...
    'pos',[hpcont+hsizecont*6.2 vpcont-vsizecont*6.4  hpcont+hsizecont/1.5 vpcont-(vpcont-vsizecont)],...
    'Tag', 'feedback' ,'Value',h.feedback,'HandleVisibility','off');

uicontrol('Style','checkbox','String','Repeat 1th trial',...
    'pos',[hpcont+hsizecont*6.2 vpcont-vsizecont*8  hpcont+hsizecont/1.5 vpcont-(vpcont-vsizecont)],...
    'Tag', 'repeatft' ,'Value',h.repeatft,'HandleVisibility','off');

uipanel('Title','Demographic data',...
    'BackgroundColor',[.9 .9 .9],...
    'Position', [.01 .56 .98 .22],...
    'FontSize',8,'FontAngle','italic');

uipanel('Title','Experiment features','FontSize',8,...
    'BackgroundColor',[.9 .9 .9],...
    'Position', [.01 .4 .98 .15],...
    'FontSize',8,'FontAngle','italic');

uipanel('Title','Hypotheses features','FontSize',8,...
    'BackgroundColor',[.9 .9 .9],...
    'Position', [.01 .26 .98 .13],...
    'FontSize',8,'FontAngle','italic');

uicontrol(mainfigure, 'Style','text',...
    'Position',[hpcont vpcont-vsizecont*11  hpcont+hsizecont/2 vpcont-(vpcont-vsizecont)],...
    'String', 'slope (beta)');

uicontrol(mainfigure, 'Style', 'edit',...
    'String', h.beta,...
    'Position',[hpcont vpcont-vsizecont*12  hpcont+hsizecont/2 vpcont-(vpcont-vsizecont)],...
    'BackgroundColor', [1 1 1],...
    'Tag', 'beta');

uicontrol(mainfigure, 'Style', 'text',...
    'Position',[hpcont+hsizecont*1.5 vpcont-vsizecont*11  hpcont+hsizecont/2 vpcont-(vpcont-vsizecont)],...
    'String', 'first midpoint');

uicontrol(mainfigure, 'Style', 'edit',...
    'Position',[hpcont+hsizecont*1.5 vpcont-vsizecont*12  hpcont+hsizecont/2 vpcont-(vpcont-vsizecont)],...
    'String', h.firstmidpoint,...
    'BackGroundColor', [1 1 1],...
    'Tag', 'firstmidpoint');

uicontrol(mainfigure, 'Style','text',...
    'Position',[hpcont+hsizecont*2.5 vpcont-vsizecont*11  hpcont+hsizecont/2 vpcont-(vpcont-vsizecont)],...
    'String', 'last midpoint');

uicontrol(mainfigure, 'Style', 'edit',...
    'String', h.lastmidpoint,...
    'Position',[hpcont+hsizecont*2.5 vpcont-vsizecont*12  hpcont+hsizecont/2 vpcont-(vpcont-vsizecont)],...
    'BackgroundColor', [1 1 1],...
    'Tag', 'lastmidpoint');

uicontrol(mainfigure, 'Style','text',...
    'Position',[hpcont+hsizecont*3.5 vpcont-vsizecont*11  hpcont+hsizecont/2 vpcont-(vpcont-vsizecont)],...
    'String', 'n. of hypotheses');

uicontrol(mainfigure, 'Style', 'edit',...
    'Position',[hpcont+hsizecont*3.5 vpcont-vsizecont*12  hpcont+hsizecont/2 vpcont-(vpcont-vsizecont)],...
    'String', h.hypnr,...
    'BackgroundColor', [1 1 1],...
    'Tag', 'hypnr');

uicontrol('Style','checkbox','String','Log scale',...
    'pos',[hpcont+hsizecont*4.5 vpcont-vsizecont*12  hpcont+hsizecont/2 vpcont-(vpcont-vsizecont)],...
    'Tag', 'logscale' ,'Value',h.islog,'HandleVisibility','off','CallBack', {@Sweet_Callback});

uicontrol(mainfigure, 'Style','text',...
    'Position',[hpcont+hsizecont*6 vpcont-vsizecont*11  hpcont+hsizecont/2 vpcont-(vpcont-vsizecont)],...
    'String', 'lambda');

uicontrol(mainfigure, 'Style', 'edit',...
    'Position',[hpcont+hsizecont*6 vpcont-vsizecont*12  hpcont+hsizecont/2 vpcont-(vpcont-vsizecont)],...
    'String', h.lambda,...
    'BackgroundColor', [1 1 1],...
    'Tag', 'lambda');

uicontrol(mainfigure, 'Style', 'text',...
    'String', 'p target',...
    'Position',[hpcont vpcont-vsizecont*14  hpcont+hsizecont/1.5 vpcont-(vpcont-vsizecont)]);

uicontrol(mainfigure, 'Style', 'edit',...
    'String', h.p_target,...
    'BackgroundColor', [1 1 1],...
    'Position',[hpcont vpcont-vsizecont*15  hpcont+hsizecont/1.5 vpcont-(vpcont-vsizecont)],...
    'Tag', 'p_target');

uicontrol('Style','checkbox','String','Sweetpoint',...
    'pos',[hpcont+hsizecont*1.2 vpcont-vsizecont*15  hpcont+50 vpcont-(vpcont-vsizecont)],...
    'Tag', 'swp' ,'Value',h.swp,'HandleVisibility','off','CallBack', {@Sweet_Callback});




buttongroup2 = uibuttongroup('visible','off');

logisticbutton = uicontrol('Style','Radio','String','logistic',...
    'pos',[hpcont+hsizecont*2.2 vpcont-vsizecont*15  hpcont+hsizecont/1.5 vpcont-(vpcont-vsizecont)],...
    'Tag', 'logistic',...
    'parent',buttongroup2,'HandleVisibility','off');

gaussianbutton = uicontrol('Style','Radio','String','gaussian',...
    'pos',[hpcont+hsizecont*2.2 vpcont-vsizecont*16  hpcont+hsizecont/1.5 vpcont-(vpcont-vsizecont)],...
    'Tag', 'gaussian',...
    'parent',buttongroup2,'HandleVisibility','off');

set(buttongroup2,'SelectionChangeFcn',@whichfunction);
set(buttongroup2,'Visible','on');






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


uiwait


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Sweet_Callback(a,b)
global h;
mh = guihandles(gcf);
logicswp = get(mh.swp,'Value') == get(mh.swp,'Max');
nafc = str2double(get(mh.nafc,'String'));
h.swp=1;
if nafc
    gamma = 1/nafc;
end

if logicswp
    if  nafc
        swp = (2*gamma+1+sqrt(1+8*gamma))/(3+sqrt(1+8*gamma));
        set(mh.p_target,'String', swp);
    else
        set(mh.p_target,'String', '0.631');
    end
else     %no sweetpoint
    if nafc
        set(mh.p_target,'String', (1+gamma)/2);
    else
        set(mh.p_target,'String', '0.5');
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function selcbk(source,eventdata)
mh = guihandles(gcf);
sel= get(eventdata.NewValue,'String');

if iscell(get(mh.swp,'Value'))
    tempVal=get(mh.swp,'Value');
    tempMax=get(mh.swp,'Max');
    logicswp= tempVal{1}==tempMax{1};

else
    logicswp=get(mh.swp,'Value')==get(mh.swp,'Max');

end

whichtask(sel,logicswp)

%%%%%%%%%%%%%%%%%%
function ChangeDefault_Callback(a,b)
global h;
savestart (1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function whichfunction(source,eventdata)

mh = guihandles(gcf);
sel= get(eventdata.NewValue,'String');
h.functiontype=sel;




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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
h.ntrials = str2double(get(mh.ntrials,'String'));
h.nblocks = str2double(get(mh.nblocks,'String'));
h.startinglevel = str2double(get(mh.startinglevel,'String'));
h.standard = str2double(get(mh.standard,'String'));
h.beta = str2double(get(mh.beta,'String'));
h.firstmidpoint = str2double(get(mh.firstmidpoint,'String'));
h.hypnr = str2double(get(mh.hypnr,'String'));
h.lastmidpoint = str2double(get(mh.lastmidpoint,'String'));
h.lambda = str2double(get(mh.lambda,'String'));
h.p_target = str2double(get(mh.p_target,'String'));
h.nafc = str2double(get(mh.nafc,'String'));
h.nafcbutton=get(mh.swp,'Value');


if get(mh.logistic,'Value')
    h.functiontype= 'logistic';
end
if get(mh.gaussian,'Value')
    h.functiontype= 'gaussian';
end





if iscell(get(mh.swp,'Value'))
    tempVal=get(mh.swp,'Value');
    tempMax=get(mh.swp,'Max');
    h.swp= tempVal{1}==tempMax{1};

else
    h.swp=get(mh.swp,'Value')==get(mh.swp,'Max');
end

if iscell(get(mh.feedback,'Value'))
    tempVal=get(mh.feedback,'Value');
    tempMax=get(mh.feedback,'Max');
    h.feedback= tempVal{1}==tempMax{1};
    tempVal2=get(mh.repeatft,'Value');
    tempMax2=get(mh.repeatft,'Max');
    h.repeatft= tempVal2{1}==tempMax2{1};
else
    h.feedback = get(mh.feedback,'Value')==get(mh.feedback,'Max');
    h.repeatft = get(mh.repeatft,'Value')==get(mh.repeatft,'Max');
end

if ~ h.nafc && h.feedback
    msgbox('No feedback can be displayed for yes/no tasks.','MLProcedure','Warn');
    set(mh.feedback,'Value', 0);
    return

end
h.description = get(mh.description,'String');

if h.firstmidpoint > h.lastmidpoint
    errordlg('Last midpoint is smaller than first midpoint','ML Procedure');

    return;
end

if iscell(get(mh.logscale,'Value'))
    tempVal=get(mh.logscale,'Value');
    tempMax=get(mh.logscale,'Max');
    h.feedback= tempVal{1}==tempMax{1};
else
    islog = get(mh.logscale,'Value') == get(mh.logscale,'Max');
end

if islog
    if h.firstmidpoint < 0 || h.lastmidpoint < 0
        sel= questdlg('Negative Logarithm! Would you like a non-log scale?',...
            'ML Procedure',...
            'Yes','No','Yes');
        switch sel,
            case 'Yes',
                h.midpoints=h.firstmidpoint:(h.lastmidpoint-h.firstmidpoint)/(h.hypnr-1):h.lastmidpoint;
            case 'No'
                return
        end
    else
        h.midpoints= 10.^(log10(h.firstmidpoint):log10(h.lastmidpoint/h.firstmidpoint)/(h.hypnr-1):log10(h.lastmidpoint));
    end
else
    h.midpoints=h.firstmidpoint:(h.lastmidpoint-h.firstmidpoint)/h.hypnr:h.lastmidpoint;
end

h.islog = islog;

if  strcmp(get(mh.textnafc,'visible'), 'on')
    h.nafc = str2double(get(mh.nafc,'String'));
    h.gamma = 1/h.nafc;
    h.catchtrial = 0;
else
    h.nafc = 0;
    h.gamma = [0, 0.1, 0.2, 0.3, 0.4];
    h.catchtrial =  str2double(get(mh.catchprop,'String'));
    h.feedback = 0;
end

if savedata
    mddefaults = ([matlabroot myslash 'toolbox' myslash 'psychoacoustics' myslash 'Defaults_MLP']);
    if exist(mddefaults, 'dir')
        defaults = [mddefaults ,myslash];
    else
        defaults = [cd ,myslash 'Defaults_MLP' myslash];
    end

    filexppar = [defaults ,h.experiment(1:end-4),'Par_MLP.m'];

    if  exist(filexppar,'file')

        fileexist=1;
        answ = questdlg('Default parameters for this experiment already exists. Do you want to replace it?', 'ML parameters','Yes','No','Yes');

    else
        fileexist=0;
    end


    if ~fileexist || strcmp(answ, 'Yes')

        logicswp=get(mh.swp,'Value')==get(mh.swp,'Max');
        fid3 =    fopen(filexppar, 'w');
        filepar = [h.experiment(1:end-4), 'Par_MLP(mh)'];
        fprintf(fid3, '%s\t','function ');
        fprintf(fid3, '%s\n',filepar);
        fprintf(fid3, '%s\t','set(mh.nblocks,''String'',');
        fprintf(fid3, '%g\t',h.nblocks);
        fprintf(fid3, '%s\n',');');
        fprintf(fid3, '%s','set(mh.fileout,''String'',''');
        fprintf(fid3, '%s',h.fileout);
        fprintf(fid3, '%s\n',''');');
        fprintf(fid3, '%s\t','set(mh.beta,''String'',');
        fprintf(fid3, '%g\t',h.beta);
        fprintf(fid3, '%s\n',');');
        fprintf(fid3, '%s\t','set(mh.ntrials,''String'',');
        fprintf(fid3, '%g\t',h.ntrials);
        fprintf(fid3, '%s\n',');');
        fprintf(fid3, '%s\t','set(mh.hypnr,''String'',');
        fprintf(fid3, '%g\t',h.hypnr);
        fprintf(fid3, '%s\n',');') ;
        fprintf(fid3, '%s\t','set(mh.firstmidpoint,''String'',');
        fprintf(fid3, '%g\t',h.firstmidpoint);
        fprintf(fid3, '%s\n',');') ;
        fprintf(fid3, '%s\t','set(mh.lastmidpoint,''String'',');
        fprintf(fid3, '%g\t',h.lastmidpoint);
        fprintf(fid3, '%s\n',');') ;
        fprintf(fid3, '%s\t','set(mh.startinglevel,''String'',');
        fprintf(fid3, '%g\t',h.startinglevel);
        fprintf(fid3, '%s\n',');') ;
        fprintf(fid3, '%s\t','set(mh.standard,''String'',');
        fprintf(fid3, '%g\t',h.standard);
        fprintf(fid3, '%s\n',');') ;
        fprintf(fid3, '%s\t','set(mh.logscale, ''value'' , ');
        fprintf(fid3, '%g\t',h.islog);
        fprintf(fid3, '%s\n',');');
        fprintf(fid3, '%s\t','set(mh.feedback, ''value'' , ');
        fprintf(fid3, '%g\t',h.feedback);
        fprintf(fid3, '%s\n',');');
        fprintf(fid3, '%s\t','set(mh.repeatft, ''value'' , ');
        fprintf(fid3, '%g\t',h.repeatft);
        fprintf(fid3, '%s\n',');');
        fprintf(fid3, '%s\t','set(mh.lambda,''String'',');
        fprintf(fid3, '%g\t',h.lambda);
        fprintf(fid3, '%s\n',');') ;
        fprintf(fid3, '%s\t','set(mh.SaveResults, ''value'' , ');
        fprintf(fid3, '%g\t',h.SaveResults);
        fprintf(fid3, '%s\n',');');
        if h.nafc
            fprintf(fid3, '%s\n','set(mh.nafcbutton, ''value'' ,1); ');
            fprintf(fid3, '%s\t','set(mh.nafc, ''String'' , ');
            fprintf(fid3, '%g\t',h.nafc);
            fprintf(fid3, '%s\n',');')  ;
            fprintf(fid3, '%s\t','whichtask(''nAFC'',');
            fprintf(fid3, '%g\t',logicswp);
            fprintf(fid3, '%s\n',');');
        else
            fprintf(fid3, '%s\t','set(mh.catchprop, ''String'' , ');
            fprintf(fid3, '%g\t',h.catchtrial);
            fprintf(fid3, '%s\n',');');
            fprintf(fid3, '%s\n','set(mh.nafcbutton, ''value'' ,0); ');
            fprintf(fid3, '%s\t','whichtask(''yes/no'',');
            fprintf(fid3, '%g\t',logicswp);
            fprintf(fid3, '%s\n',');');
            fprintf(fid3, '%s\t','set(mh.yes_nobutton, ''value'' ,1); ');
            fprintf(fid3, '%s\t','set(mh.nafc, ''String'' , 0);')  ;
        end

        fclose(fid3);
    end

    return
else

    h.StimulusLevel= zeros(1,h.ntrials );
    h.SubjectAccuracy= zeros(1,h.ntrials );
    h.FA = zeros(1,h.ntrials );
    h.TemporaryThreshold=zeros(1,h.ntrials );
    h.MATSAVEDATA = zeros(h.nblocks*h.ntrials,6);
    h.StimulusLevel(1)=h.startinglevel;   %first trial: stimuli intensity = assigned intensity
end;

delete (gcf);
if gcf
    delete (gcf);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function closeGUImasolo(src,evnt)
global h;

selection = questdlg('Do you want to close MLParameters?',...
    'ML Procedure',...
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Updateptarget_Callback(source,eventdata)
global h;
mh = guihandles(gcf);
h.nafc = str2double(get(mh.nafc,'String'));
logicswp = get(mh.swp,'Value') == get(mh.swp,'Max');
gamma = 1/h.nafc;
if gamma == Inf
    gamma = 0;
end

if logicswp
    swp = (2*gamma+1+sqrt(1+8*gamma))/(3+sqrt(1+8*gamma));
    set(mh.p_target,'String', swp);
else
    set(mh.p_target,'String', (1+gamma)/2);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function whichtask (sel,logicswp)
global h
mh = guihandles(gcf);
switch sel
    case 'yes/no'
        h.tasktype = 0;
        set(mh.textcatchtrial,'Visible', 'on');
        set(mh.catchprop,'Visible', 'on');
        set(mh.nafc,'Visible', 'off');
        set(mh.textnafc,'Visible','off');
        set(mh.nafc,'String', 0);
        set(mh.feedback,'Value', 0);
        if logicswp
            set(mh.p_target,'String', '0.631');
            set (mh.swp,'Value',1);
        else
            set(mh.p_target,'String', '0.5');
            set (mh.swp,'Value',0);
        end
    case 'nAFC'
        h.tasktype = 1;
        nafc= str2double(get(mh.nafc,'String'));
        if nafc
            set(mh.nafc,'String',nafc);
        else
            set(mh.nafc,'String',2);
        end
        set(mh.ptargetupdate,'Visible', 'on');
        gamma = 1/str2double( get(mh.nafc,'String'));
        set(mh.nafc,'Visible', 'on');
        set(mh.textnafc,'Visible','on');
        set(mh.textcatchtrial,'Visible', 'of');
        set(mh.catchprop,'Visible', 'of');

        if logicswp
            swp = (2*gamma+1+sqrt(1+8*gamma))/(3+sqrt(1+8*gamma));
            set(mh.p_target,'String', swp);
            set (mh.swp,'Value',1);
        else
            set(mh.p_target,'String', (1+gamma)/2);
            set (mh.swp,'Value',0);
        end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function exp(source,eventdata,indx)
global h;
global d;
global myslash;
mddefaults = ([matlabroot myslash 'toolbox' myslash 'psychoacoustics' myslash 'Defaults_MLP']);

if exist(mddefaults, 'dir')
    defaults = [mddefaults ,myslash];
else
    defaults = [cd ,myslash 'Defaults_MLP' myslash];
end

filexp = [defaults ,'lastexp.txt'];

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
        warndlg('You are running the procedure outside the psychoacustics directory. Either enter the parameters or change directory and try again','MLProcedure procedure');
    end
end
mh = guihandles(gcf);


if size(d,1)< indx
    indx=1;
    msgbox('Last experiment defaults have probably been deleted. The defaults refer to the first experiment in the list.','MLProcedure','Warn')
end;
filesel = d(indx).name ;


h.experiment = filesel(1:length(filesel)-2) ;
set(mh.ExpLabel,'String', h.experiment);
descr = help(h.experiment);
set(mh.description,'String', descr);
funpar = [h.experiment(1:end-4), 'Par_MLP(mh)'];
filepar= [h.experiment(1:end-4), 'Par_MLP.m'];

if  exist(filepar,'file')
    eval(funpar);
else
    set(mh.ntrials,'String','');
    set(mh.nblocks,'String','');
    set(mh.startinglevel,'String','');
    set(mh.standard,'String','');
    set(mh.nafc,'Value',0);
    set(mh.p_target ,'String','');
    set(mh.beta ,'String','');
    set(mh.firstmidpoint  ,'String','');
    set(mh.lastmidpoint   ,'String','');
    set(mh.hypnr   ,'String','');
    set(mh.feedback   ,'Value',0);
    set(mh.lambda   ,'String','');
    set(mh.fileout   ,'String','');
    set(mh.SaveResults   ,'Value',0);
    set(mh.swp   ,'Value',0);
    set(mh.logscale   ,'Value',0);
    warndlg(['File ' filepar ' doesn''t exist. Press SAVE DAFAULTS if you want to save defaults for this experiment'], 'ML parameters');
end
fid2 =    fopen(filexp, 'w');
fprintf(fid2, '%g', indx);
fclose(fid2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function editexp(source,eventdata,indx)
global d;
filesel = d(indx).name;
edit (filesel);

