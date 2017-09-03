function 	ITD_DiscriminationPureTonePar_MLP(mh)
set(mh.nblocks,'String',	5	);
set(mh.fileout,'String','datafile.txt');
set(mh.beta,'String',	100	);
set(mh.ntrials,'String',	30	);
set(mh.hypnr,'String',	100	);
set(mh.firstmidpoint,'String',	0.0001	);
set(mh.lastmidpoint,'String',	0.3	);
set(mh.startinglevel,'String',	0.3	);
set(mh.standard,'String',	0	);
set(mh.logscale, 'value' , 	1	);
set(mh.feedback, 'value' , 	1	);
set(mh.repeatft, 'value' , 	1	);
set(mh.lambda,'String',	0	);
set(mh.SaveResults, 'value' , 	1	);
set(mh.nafcbutton, 'value' ,1); 
set(mh.nafc, 'String' , 	2	);
whichtask('nAFC',	1	);
