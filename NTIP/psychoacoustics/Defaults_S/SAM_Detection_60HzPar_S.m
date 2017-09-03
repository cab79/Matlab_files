function 	SAM_Detection_60HzPar_S(mh)
set(mh.nblocks,'String',	5	);
set(mh.fileout,'String','data.txt');
set(mh.startinglevel,'String',	0	);
set(mh.standard,'String',	-Inf	);
set(mh.reversals,'String',	'4 8 ');
set(mh.stepsize,'String',	'3 1.5 ');
set(mh.NameStepsize,'String', 'Step Size');
set(mh.stepbutton,'value',1	);
set(mh.TwoDownOneUp,'value',1	);
set(mh.Arithmetic,'value',1	);
set(mh.nafcbutton, 'value' ,1); 
set(mh.textnafc,'Visible','on'); 
set(mh.nafc,'Visible', 'on'); 
set(mh.nafc, 'String' , 	3	);
set(mh.reversalForthresh,'String',	8	);
set(mh.feedback, 'value' , 	1	);
set(mh.SaveResults, 'value' , 	1	);
