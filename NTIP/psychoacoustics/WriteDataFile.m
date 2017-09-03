function WriteDataFile(tracking, myfile, MATSAVEDATA, nsub, name, gender, age, note, thresholdsonly)
%WRITEDATAFILE SAVES THE DATA INTO TWO FILES.THE FIRST, WHICH NAME IS GIVEN BY
%THE EXPERIMENTER CONTAINS ALL THE RAW DATA, THE SECOND, THRESHOLD.TXT,
%CONTAINS THE THRESHOLDS FOR EACH SUBJECT FOR EACH BLOCK

if nargin<9
    thresholdsonly = 0;
end

%STRUCT = ML PROCEDURE STRUCTURE
if isempty (myfile)
    myfile = 'data.txt';
end
% fileout = ([cd, '\', myfile]);
fileout = myfile;

if isempty(name)
    namesubj = 'untitled.txt';
else
    namesubj = [name '.txt'];
end
filesubj = fopen(namesubj, 'a');
fprintf (filesubj, 'BLOCK\tTHRESHOLD\n');
if exist(fileout, 'file');
    datafile = fopen(fileout, 'a');
else
    datafile = fopen(fileout, 'a');
    switch tracking
        case 'MLP'
            fprintf (datafile, 'nsub\tname\tsex\tage\tblock\ttrial\tlevel\tgamma\trispAC\tthreshold\tsubnote\n');
        case 'Staircase'
            fprintf (datafile, 'nsub\tname\tsex\tage\tblock\ttrial\tlevel\trispAC\tnreversals\tactualstep\tsubnote\n');
        case 'PEST'
            fprintf (datafile, 'nsub\tname\tsex\tage\tblock\ttrial\tlevel\trispAC\tstep\tsubnote\n');
    end
end
switch tracking
    case 'MLP'
        for i=1:size(MATSAVEDATA,1)
            fprintf(datafile, '%3.0f\t%s\t%s\t%2.0f\t%2.0f\t%2.0f\t%4.3f\t%1.1f\t%1.0f\t%4.3f\t%s\n',nsub,name,gender,age,...
                MATSAVEDATA(i,1),MATSAVEDATA(i,2),MATSAVEDATA(i,3),MATSAVEDATA(i,4),MATSAVEDATA(i,5),...
                MATSAVEDATA(i,6),note);
            if i<size(MATSAVEDATA,1)&& MATSAVEDATA (i+1,1)~=  MATSAVEDATA (i,1)|| i==size(MATSAVEDATA,1)
                fprintf (filesubj,'%2.0f\t%3.3f\n',MATSAVEDATA (i,1),MATSAVEDATA (i,6));
            end
        end
    case 'Staircase'
        for i=1:size(MATSAVEDATA,1)
            fprintf(datafile, '%3.0f\t%s\t%s\t%2.0f\t%2.0f\t%2.0f\t%4.3f\t%1.1f\t%1.0f\t%4.3f\t%s\n',nsub,name,gender,age,...
                MATSAVEDATA(i,1),MATSAVEDATA(i,2),MATSAVEDATA(i,3),MATSAVEDATA(i,4),MATSAVEDATA(i,5),...
                MATSAVEDATA(i,6),note);
            if i<size(MATSAVEDATA,1)&& MATSAVEDATA (i+1,1)~=  MATSAVEDATA (i,1)|| i==size(MATSAVEDATA,1)
                fprintf (filesubj,'%2.0f\t%3.3f\n',MATSAVEDATA (i,1),thresholdsonly(MATSAVEDATA (i,1)));
            end
        end
    case 'PEST'
        for i=1:size(MATSAVEDATA,1)
            fprintf(datafile, '%3.0f\t%s\t%s\t%2.0f\t%2.0f\t%2.0f\t%4.3f\t%1.0f\t%4.3f\t%s\n',nsub,name,gender,age,...
                MATSAVEDATA(i,1),MATSAVEDATA(i,2),MATSAVEDATA(i,3),MATSAVEDATA(i,4),MATSAVEDATA(i,5),note);
            if i<size(MATSAVEDATA,1)&& MATSAVEDATA (i+1,1)~=  MATSAVEDATA (i,1)|| i==size(MATSAVEDATA,1)
                fprintf (filesubj,'%2.0f\t%3.3f\n',MATSAVEDATA (i,1),thresholdsonly(MATSAVEDATA (i,1)));
            end
        end
end

fclose (datafile);
fclose (filesubj);
if  strcmp (namesubj,'untitled.txt')
    msgbox('Subject''s name missed. Thresholds have been saved in the "untitled.txt" file.','Psychoacoustics','Warn')
end

if  strcmp ( myfile, 'data.txt')
    msgbox('File data name is missing. Data has been saved in the "data.txt" file.','Psychoacoustics','Warn')
end
