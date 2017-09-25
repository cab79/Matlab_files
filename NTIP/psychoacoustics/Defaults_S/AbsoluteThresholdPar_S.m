function 	AbsoluteThresholdPar_S(mh)
set(mh.nblocks,'String',	5	);
set(mh.fileout,'String','data.txt');
set(mh.startinglevel,'String',	-50	);
set(mh.standard,'String',	-Inf	);
set(mh.reversals,'String',	'4 8 ');
set(mh.stepsize,'String',	'2 1 ');
set(mh.NameStepsize,'String', 'Step Size');
set(mh.stepbutton,'value',1	);
set(mh.SimpleUpdown,'value',1	);
set(mh.Arithmetic,'value',1	);
set(mh.yes_nobutton, 'value' ,1); 
set(mh.nafc, 'String' , 0);
set(mh.reversalForthresh,'String',	8	);
set(mh.feedback, 'value' , 	0	);
set(mh.SaveResults, 'value' , 	1	);
