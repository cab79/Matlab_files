function 	PitchDiscriminationPureTonePar_S(mh)
set(mh.nblocks,'String',	5	);
set(mh.fileout,'String','data.txt');
set(mh.startinglevel,'String',	100	);
set(mh.standard,'String',	1000	);
set(mh.reversals,'String',	'4 8 ');
set(mh.stepsize,'String',	'2 1.41421 ');
set(mh.NameStepsize,'String', '	Factor ');
set(mh.factorbutton,'value',1	);
set(mh.TwoDownOneUp,'value',1	);
set(mh.Arithmetic,'value',1	);
set(mh.nafcbutton, 'value' ,1); 
set(mh.textnafc,'Visible','on'); 
set(mh.nafc,'Visible', 'on'); 
set(mh.nafc, 'String' , 	3	);
set(mh.reversalForthresh,'String',	8	);
set(mh.feedback, 'value' , 	1	);
set(mh.SaveResults, 'value' , 	1	);
