function 	ITD_DiscriminationSAM_PureTonePar_Pest(mh)
set(mh.nblocks,'String',	5	);
set(mh.fileout,'String','data.txt');
set(mh.startinglevel,'String',	0.7	);
set(mh.standard,'String',	NaN	);
set(mh.minlevel,'String',	0	);
set(mh.nafc, 'String' , 	3	);
set(mh.feedback, 'value' , 	1	);
set(mh.p_target,'String',	0.75	);
set(mh.waldfactor,'String',	1	);
set(mh.startingstepsize,'String',	0.2	);
set(mh.finalstepsize,'String',	0.025	);
set(mh.maxstepsize,'String',	0.4	);
set(mh.SaveResults, 'value' , 	1	);
