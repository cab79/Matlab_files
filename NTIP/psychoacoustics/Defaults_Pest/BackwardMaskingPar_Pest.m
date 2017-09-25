function 	BackwardMaskingPar_Pest(mh)
set(mh.nblocks,'String',	5	);
set(mh.fileout,'String','data.txt');
set(mh.startinglevel,'String',	-20	);
set(mh.standard,'String',	-Inf	);
set(mh.minlevel,'String',	-Inf	);
set(mh.nafc, 'String' , 	3	);
set(mh.feedback, 'value' , 	1	);
set(mh.p_target,'String',	0.75	);
set(mh.waldfactor,'String',	1	);
set(mh.startingstepsize,'String',	8	);
set(mh.finalstepsize,'String',	0.5	);
set(mh.maxstepsize,'String',	16	);
set(mh.SaveResults, 'value' , 	1	);
