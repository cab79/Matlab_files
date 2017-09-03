function 	TemporalOrderTonesPar_MLP(mh)
set(mh.nblocks,'String',	5	);
set(mh.fileout,'String','data.txt');
set(mh.beta,'String',	0.1	);
set(mh.ntrials,'String',	30	);
set(mh.hypnr,'String',	100	);
set(mh.firstmidpoint,'String',	20	);
set(mh.lastmidpoint,'String',	300	);
set(mh.startinglevel,'String',	200	);
set(mh.standard,'String',	NaN	);
set(mh.logscale, 'value' , 	1	);
set(mh.feedback, 'value' , 	0	);
set(mh.repeatft, 'value' , 	1	);
set(mh.lambda,'String',	0	);
set(mh.SaveResults, 'value' , 	1	);
set(mh.nafcbutton, 'value' ,1); 
set(mh.nafc, 'String' , 	2	);
whichtask('nAFC',	1	);
