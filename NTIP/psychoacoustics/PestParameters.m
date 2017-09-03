function h = PestParameters()

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
h.standard = '';
h.minlevel = '';
h.nafc = '';

h.feedback= 0;


h.fileout = '';
h.SaveResults =1;
h.tasktype = 0;


h.exppos = 1;
h.experiment = '';
h.description = '';

h.p_target = '';
h.waldfactor = '';


h.finalstepsize = '';
h.startingstepsize = '';
h.maxstepsize = '';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
mainfigure = figure (100);
set(0,'Units','pixels')
scnsize = get(0, 'ScreenSize');
hpcont = scnsize(3)/36;                     %horizontal posizion control
hsizecont = scnsize(3)/14.2;                %horizontal size control
vpcont = scnsize(4)/2;                      %vertical posizion control
vsizecont = scnsize(4)/40;                  %vertical size control
md = ([matlabroot myslash 'toolbox' myslash 'psychoacoustics' myslash 'Experiments_Pest']);

if exist(md, 'dir')
    d = dir([md myslash '*.m']);
else
    d = dir([cd myslash 'Experiments_Pest'  myslash '*.m']);
end;
if ~exist('d','var')
    errordlg('Directory Experiments not found','Pest Procedure')
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
    'Name', ' Pest Procedure');

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
    'String', 'n. of blocks');

uicontrol(mainfigure, 'Style', 'edit',...
    'String', h.nblocks,...
    'Position',[hpcont vpcont-vsizecont*8  hpcont+hsizecont/2 vpcont-(vpcont-vsizecont)],...
    'BackgroundColor', [1 1 1],...
    'Tag', 'nblocks');

uicontrol(mainfigure, 'Style', 'text',...
    'String', 'starting level',...
    'Position',[hpcont+hsizecont*1.2 vpcont-vsizecont*7  hpcont+hsizecont/2 vpcont-(vpcont-vsizecont)]);

uicontrol(mainfigure, 'Style', 'edit',...
    'String', h.startinglevel,...
    'BackgroundColor', [1 1 1],...
    'Position',[hpcont+hsizecont*1.2 vpcont-vsizecont*8  hpcont+hsizecont/2 vpcont-(vpcont-vsizecont)],...
    'Tag', 'startinglevel');

uicontrol(mainfigure, 'Style', 'text',...
    'String', 'standard level',...
    'Position',[hpcont+hsizecont*2.4 vpcont-vsizecont*7  hpcont+hsizecont/2 vpcont-(vpcont-vsizecont)]);

uicontrol(mainfigure, 'Style', 'edit',...
    'String', h.standard,...
    'BackgroundColor', [1 1 1],...
    'Position',[hpcont+hsizecont*2.4 vpcont-vsizecont*8  hpcont+hsizecont/2 vpcont-(vpcont-vsizecont)],...
    'Tag', 'standard');

uicontrol(mainfigure, 'Style','text',...
    'String', 'minimum level',...
    'Position',[hpcont+hsizecont*3.6 vpcont-vsizecont*7  hpcont+hsizecont/2 vpcont-(vpcont-vsizecont)]);

uicontrol(mainfigure, 'Style', 'edit',...
    'String', h.minlevel,...
    'Position',[hpcont+hsizecont*3.6 vpcont-vsizecont*8  hpcont+hsizecont/2 vpcont-(vpcont-vsizecont)],...
    'BackgroundColor', [1 1 1],...
    'Tag', 'minlevel');

uicontrol(mainfigure, 'Style','text',...
    'String', 'nAFC',...
    'Position',[hpcont+hsizecont*4.8 vpcont-vsizecont*7  hpcont+hsizecont/2 vpcont-(vpcont-vsizecont)]);

uicontrol(mainfigure, 'Style', 'edit',...
    'String', h.nafc,...
    'Position',[hpcont+hsizecont*4.8 vpcont-vsizecont*8  hpcont+hsizecont/2 vpcont-(vpcont-vsizecont)],...
    'BackgroundColor', [1 1 1],...
    'Tag', 'nafc');

uicontrol('Style','checkbox','String','Feedback',...
    'pos',[hpcont+hsizecont*6.0 vpcont-vsizecont*8  hpcont+hsizecont/2 vpcont-(vpcont-vsizecont)],...
    'Tag', 'feedback' ,'Value',h.feedback,'HandleVisibility','off');

uicontrol(mainfigure, 'Style','text',...
    'Position',[hpcont vpcont-vsizecont*11  hpcont+hsizecont/2 vpcont-(vpcont-vsizecont)],...
    'String', 'p_target');

uicontrol(mainfigure, 'Style', 'edit',...
    'String', h.p_target,...
    'Position',[hpcont vpcont-vsizecont*12  hpcont+hsizecont/2 vpcont-(vpcont-vsizecont)],...
    'BackgroundColor', [1 1 1],...
    'Tag', 'p_target');

uicontrol(mainfigure, 'Style', 'text',...
    'Position',[hpcont+hsizecont*1.2 vpcont-vsizecont*11  hpcont+hsizecont/2 vpcont-(vpcont-vsizecont)],...
    'String', 'W constant');

uicontrol(mainfigure, 'Style', 'edit',...
    'Position',[hpcont+hsizecont*1.2 vpcont-vsizecont*12  hpcont+hsizecont/2 vpcont-(vpcont-vsizecont)],...
    'String', h.waldfactor,...
    'BackGroundColor', [1 1 1],...
    'Tag', 'waldfactor');

uicontrol(mainfigure, 'Style','text',...
    'Position',[hpcont+hsizecont*2.4 vpcont-vsizecont*11  hpcont+hsizecont/2 vpcont-(vpcont-vsizecont)],...
    'String', 'Starting step size');


uicontrol(mainfigure, 'Style', 'edit',...
    'String', h.startingstepsize,...
    'Position',[hpcont+hsizecont*2.4 vpcont-vsizecont*12  hpcont+hsizecont/2 vpcont-(vpcont-vsizecont)],...
    'BackgroundColor', [1 1 1],...
    'Tag', 'startingstepsize');



uicontrol(mainfigure, 'Style','text',...
    'Position',[hpcont+hsizecont*3.6 vpcont-vsizecont*11  hpcont+hsizecont/2 vpcont-(vpcont-vsizecont)],...
    'String', 'Final step size');


uicontrol(mainfigure, 'Style', 'edit',...
    'String', h.finalstepsize,...
    'Position',[hpcont+hsizecont*3.6 vpcont-vsizecont*12  hpcont+hsizecont/2 vpcont-(vpcont-vsizecont)],...
    'BackgroundColor', [1 1 1],...
    'Tag', 'finalstepsize');

uicontrol(mainfigure, 'Style','text',...
    'Position',[hpcont+hsizecont*4.8 vpcont-vsizecont*11  hpcont+hsizecont/2 vpcont-(vpcont-vsizecont)],...
    'String', 'Max step size');


uicontrol(mainfigure, 'Style', 'edit',...
    'String', h.maxstepsize,...
    'Position',[hpcont+hsizecont*4.8 vpcont-vsizecont*12  hpcont+hsizecont/2 vpcont-(vpcont-vsizecont)],...
    'BackgroundColor', [1 1 1],...
    'Tag', 'maxstepsize');


uipanel('Title','Demographic data',...
    'BackgroundColor',[.9 .9 .9],...
    'Position', [.01 .56 .98 .22],...
    'FontSize',8,'FontAngle','italic');

uipanel('Title','Experiment features','FontSize',8,...
    'BackgroundColor',[.9 .9 .9],...
    'Position', [.01 .4 .98 .15],...
    'FontSize',8,'FontAngle','italic');

uipanel('Title','Pest features','FontSize',8,...
    'BackgroundColor',[.9 .9 .9],...
    'Position', [.01 .26 .98 .13],...
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


uiwait


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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
h.minlevel = str2double(get(mh.minlevel,'String'));
h.nafc = str2double(get(mh.nafc,'String'));


h.p_target = str2double(get(mh.p_target,'String'));
h.waldfactor =str2double(get(mh.waldfactor,'String'));
h.finalstepsize = str2double(get(mh.finalstepsize,'String'));
h.startingstepsize = str2double(get(mh.startingstepsize,'String'));
h.maxstepsize = str2double(get(mh.maxstepsize,'String'));


if h.startingstepsize <= h.finalstepsize
    warndlg('Starting step size cannot be smaller than final step size, please try again','Pest procedure');
    return;
end

if h.startingstepsize >= h.maxstepsize
    warndlg('Starting step size cannot be larger than max step size, please try again','Pest procedure');
    return;
end


if iscell(get(mh.feedback,'Value'))
    tempVal=get(mh.feedback,'Value');
    tempMax=get(mh.feedback,'Max');
    h.feedback= tempVal{1}==tempMax{1};

else
    h.feedback = get(mh.feedback,'Value')==get(mh.feedback,'Max');
end

h.description = get(mh.description,'String');


if savedata
    mddefaults = ([matlabroot myslash 'toolbox' myslash 'psychoacoustics' myslash 'Defaults_Pest']);
    if exist(mddefaults, 'dir')
        defaults = [mddefaults ,myslash];
    else
        defaults = [cd ,myslash 'Defaults_Pest', myslash];
    end

    filexppar = [defaults ,h.experiment(1:end-5),'Par_Pest.m'];



    if  exist(filexppar,'file')

        fileexist=1;
        answ = questdlg('Default parameters for this experiment already exists. Do you want to replace it?', 'Pest Procedure','Yes','No','Yes');

    else
        fileexist=0;
    end


    if ~fileexist || strcmp(answ, 'Yes')


        fid3 =    fopen(filexppar, 'w');
        filepar = [h.experiment(1:end-5), 'Par_Pest(mh)'];
        fprintf(fid3, '%s\t','function ');
        fprintf(fid3, '%s\n',filepar);
        fprintf(fid3, '%s\t','set(mh.nblocks,''String'',');
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

        fprintf(fid3, '%s\t','set(mh.minlevel,''String'',');
        fprintf(fid3, '%g\t',h.minlevel);
        fprintf(fid3, '%s\n',');') ;

        fprintf(fid3, '%s\t','set(mh.nafc, ''String'' , ');
        fprintf(fid3, '%g\t',h.nafc);
        fprintf(fid3, '%s\n',');')  ;

        fprintf(fid3, '%s\t','set(mh.feedback, ''value'' , ');
        fprintf(fid3, '%g\t',h.feedback);
        fprintf(fid3, '%s\n',');');


        fprintf(fid3, '%s\t','set(mh.p_target,''String'',');
        fprintf(fid3, '%g\t',h.p_target);
        fprintf(fid3, '%s\n',');') ;

        fprintf(fid3, '%s\t','set(mh.waldfactor,''String'',');
        fprintf(fid3, '%g\t',h.waldfactor);
        fprintf(fid3, '%s\n',');') ;

        fprintf(fid3, '%s\t','set(mh.startingstepsize,''String'',');
        fprintf(fid3, '%g\t',h.startingstepsize);
        fprintf(fid3, '%s\n',');') ;

        fprintf(fid3, '%s\t','set(mh.finalstepsize,''String'',');
        fprintf(fid3, '%g\t',h.finalstepsize);
        fprintf(fid3, '%s\n',');') ;

        fprintf(fid3, '%s\t','set(mh.maxstepsize,''String'',');
        fprintf(fid3, '%g\t',h.maxstepsize);
        fprintf(fid3, '%s\n',');') ;

        fprintf(fid3, '%s\t','set(mh.SaveResults, ''value'' , ');
        fprintf(fid3, '%g\t',h.SaveResults);
        fprintf(fid3, '%s\n',');');

        fclose(fid3);
    end

    return
else

    h.StimulusLevel= zeros(1,h.nblocks );
    h.SubjectAccuracy= zeros(1,h.nblocks );
    h.FA = zeros(1,h.nblocks );
    h.TemporaryThreshold=zeros(1,h.nblocks );
    h.MATSAVEDATA = zeros(h.nblocks,6);
    h.StimulusLevel(1)=h.startinglevel;   %first trial: stimuli intensity = assigned intensity
end;

delete (gcf);
if gcf
    delete (gcf);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function closeGUImasolo(src,evnt)
global h;

selection = questdlg('Do you want to close Pest Procedure?',...
    'Pest Procedure',...
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function exp(source,eventdata,indx)
global h;
global d;
global myslash;
mddefaults = ([matlabroot myslash 'toolbox' myslash 'psychoacoustics' myslash 'Defaults_Pest']);

if exist(mddefaults, 'dir')
    defaults = [mddefaults ,myslash];
else
    defaults = [cd ,myslash 'Defaults_Pest', myslash];
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

        warndlg('You are running the procedure outside the psychoacustics directory. Either enter the parameters or change directory and try again','Pest procedure');

    end
end

mh = guihandles(gcf);


if size(d,1)< indx
    indx=1;
    warndlg('Last experiment defaults have probably been deleted. The defaults refer to the first experiment in the list.','Pestrocedure','Warn')
end;
filesel = d(indx).name ;


h.experiment = filesel(1:length(filesel)-2) ;
set(mh.ExpLabel,'String', h.experiment);
descr = help(h.experiment);
set(mh.description,'String', descr);
funpar = [h.experiment(1:end-5), 'Par_Pest(mh)'];
filepar= [h.experiment(1:end-5), 'Par_Pest.m'];


if  exist(filepar,'file')
    eval(funpar);
else
    set(mh.nblocks,'String','');
    set(mh.startinglevel,'String','');
    set(mh.standard,'String','');
    set(mh.minlevel,'String','');
    set(mh.nafc,'String','');

    set(mh.feedback   ,'Value',0);

    set(mh.fileout   ,'String','');
    set(mh.SaveResults   ,'Value',0);

    warndlg(['File ' filepar ' doesn''t exist. Press SAVE DAFAULTS if you want to save defaults for this experiment'], 'Pest parameters');
end
fid2 =    fopen(filexp, 'w');
fprintf(fid2, '%g', indx);
fclose(fid2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function editexp(source,eventdata,indx)
global d;
filesel = d(indx).name;
edit (filesel);

