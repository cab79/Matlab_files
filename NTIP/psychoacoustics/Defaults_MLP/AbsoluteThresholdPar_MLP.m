function 	AbsoluteThresholdPar_MLP(mh)
set(mh.nblocks,'String',	5	);
set(mh.fileout,'String','data.txt');
set(mh.beta,'String',	1	);
set(mh.ntrials,'String',	30	);
set(mh.hypnr,'String',	100	);
set(mh.firstmidpoint,'String',	-90	);
set(mh.lastmidpoint,'String',	-10	);
set(mh.startinglevel,'String',	-10	);
set(mh.standard,'String',	-Inf	);
set(mh.logscale, 'value' , 	0	);
set(mh.feedback, 'value' , 	0	);
set(mh.repeatft, 'value' , 	1	);
set(mh.lambda,'String',	0	);
set(mh.SaveResults, 'value' , 	1	);
set(mh.catchprop, 'String' , 	0	);
set(mh.nafcbutton, 'value' ,0); 
whichtask('yes/no',	1	);
set(mh.yes_nobutton, 'value' ,1); 	set(mh.nafc, 'String' , 0);	