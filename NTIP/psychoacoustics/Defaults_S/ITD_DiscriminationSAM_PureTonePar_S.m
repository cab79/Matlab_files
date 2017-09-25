function 	ITD_DiscriminationSAM_PureTonePar_S(mh)
set(mh.nblocks,'String',	5	);
set(mh.fileout,'String','data.txt');
set(mh.startinglevel,'String',	0.8	);
set(mh.standard,'String',	NaN	);
set(mh.reversals,'String',	'4 8 ');
set(mh.stepsize,'String',	'2 1.41421 ');
set(mh.NameStepsize,'String', '');
set(mh.factorbutton,'value',1	);
set(mh.ThreeDownOneUp,'value',1	);
set(mh.Arithmetic,'value',1	);
set(mh.nafcbutton, 'value' ,1); 
set(mh.textnafc,'Visible','on'); 
set(mh.nafc,'Visible', 'on'); 
set(mh.nafc, 'String' , 	2	);
set(mh.reversalForthresh,'String',	8	);
set(mh.feedback, 'value' , 	1	);
set(mh.SaveResults, 'value' , 	1	);
