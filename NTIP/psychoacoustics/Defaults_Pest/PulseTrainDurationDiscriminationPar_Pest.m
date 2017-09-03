function 	PulseTrainDurationDiscriminationPar_Pest(mh)
set(mh.nblocks,'String',	5	);
set(mh.fileout,'String','data.txt');
set(mh.startinglevel,'String',	60	);
set(mh.standard,'String',	NaN	);
set(mh.minlevel,'String',	0	);
set(mh.nafc, 'String' , 	3	);
set(mh.feedback, 'value' , 	1	);
set(mh.p_target,'String',	0.75	);
set(mh.waldfactor,'String',	1	);
set(mh.startingstepsize,'String',	15	);
set(mh.finalstepsize,'String',	3	);
set(mh.maxstepsize,'String',	30	);
set(mh.SaveResults, 'value' , 	1	);
