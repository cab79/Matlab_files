function 	DurationDiscriminationComplexTonePar_MLP(mh)
set(mh.nblocks,'String',	5	);
set(mh.fileout,'String','data.txt');
set(mh.beta,'String',	0.1	);
set(mh.ntrials,'String',	30	);
set(mh.hypnr,'String',	100	);
set(mh.firstmidpoint,'String',	250.1	);
set(mh.lastmidpoint,'String',	450.1	);
set(mh.startinglevel,'String',	450	);
set(mh.standard,'String',	250.1	);
set(mh.logscale, 'value' , 	1	);
set(mh.feedback, 'value' , 	1	);
set(mh.repeatft, 'value' , 	1	);
set(mh.lambda,'String',	0	);
set(mh.SaveResults, 'value' , 	1	);
set(mh.nafcbutton, 'value' ,1); 
set(mh.nafc, 'String' , 	2	);
whichtask('nAFC',	1	);
