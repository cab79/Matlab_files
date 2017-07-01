clear all
subjects = {'H1', 'H2', 'H3', 'H4', 'H5', 'H6', 'H7', 'H8', 'H9', 'H10', 'H11', 'H12', 'H13' 'H14', 'H15', 'H16', 'F1', 'F2', 'F3', 'F4', 'F5', 'F6', 'F7', 'F8', 'F9', 'F10', 'F11', 'F12', 'F13', 'F14', 'F15', 'F16', 'OA1', 'OA4', 'OA5', 'OA6', 'OA7', 'OA8', 'OA9', 'OA10', 'OA11', 'OA12', 'OA13'};

fnames={'spm_epoch_';
  };

Ns = length(subjects);

for i = 1:Ns
    
subject = subjects(i);
subject = char(subject);   

fname=char(fnames);
fname= [fname subject '.mat'];
    
    D = spm_eeg_averageCAB(fullfile(pwd,fname));
    
    eval(['save ' 'fname subject '.mat' ' D']) 
    
end