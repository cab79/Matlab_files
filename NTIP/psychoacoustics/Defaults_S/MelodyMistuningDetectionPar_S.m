function 	MelodyMistuningDetectionPar_S(mh)
set(mh.nblocks,'String',	5	);
set(mh.fileout,'String','data.txt');
set(mh.startinglevel,'String',	80	);
set(mh.standard,'String',	700	);
set(mh.reversals,'String',	'4 8 ');
set(mh.stepsize,'String',	'2 1.41421 ');
set(mh.NameStepsize,'String', 'Factor');
set(mh.factorbutton,'value',1	);
set(mh.TwoDownOneUp,'value',1	);
set(mh.Arithmetic,'value',1	);
set(mh.yes_nobutton, 'value' ,1); 
set(mh.nafc, 'String' , 0);
set(mh.reversalForthresh,'String',	8	);
set(mh.feedback, 'value' , 	0	);
set(mh.SaveResults, 'value' , 	1	);
