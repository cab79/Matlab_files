function 	ProfileAnalysisPar_MLP(mh)
set(mh.nblocks,'String',	5	);
set(mh.fileout,'String','data.txt');
set(mh.beta,'String',	1	);
set(mh.ntrials,'String',	30	);
set(mh.startinglevel,'String',	-19.01	);
set(mh.standard,'String',	-40	);
set(mh.logscale, 'value' , 	0	);
set(mh.feedback, 'value' , 	1	);
set(mh.SaveResults, 'value' , 	1	);
set(mh.nafc, 'String' , 	3	);
whichtask('nAFC',	1	);
