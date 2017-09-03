function exptype = psychoacoustics()
answ=questdlg('Choose your procedure',...
    'Procedure',...
    'MLP','Staircase','Pest','MLP');
switch answ
case 'MLP'
    MLP
    case 'Staircase'
        Staircase
    case 'Pest'
        PEST
   
end