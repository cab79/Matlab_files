clear all
close all


A={
'P4_Exp2.Right.aff.flip.set'
'P7_Exp2.Left.aff.set'
'P12_Exp2.Left.aff.set'
'P15_Exp2_Right_aff.flip.set'
'P17_Exp2.Left.aff.set'
'P19_Exp2.Left.aff.set'
'P8_Exp2.Left.aff.set'
'P10_Exp2.right.aff.flip.set'
'P13_Exp2_Left.aff.set'
'P14_Exp2.right.aff.flip.set'
'P16_Exp2.Left.aff.set'
'P18_Exp2.Left.aff.set'
'P20_Exp2.Left.aff.set'
'P4_Exp2.Left.Unaff.flip.flip.set'
'P7_Exp2.Right.Unaff.flip.set'
'P12_Exp2.Right.Unaff.flip.set'
'P15_Exp2_Left_unaff.flip.flip.set'
'P17_Exp2.Right.unaff.flip.set'
'P19_Exp2.Right.Unaff.flip.set'
'P8_Exp2.Right.Unaff.flip.set'
'P10_Exp2.left.unaff.flip.flip.set'
'P13_Exp2_Right.unaff.flip.set'
'P14_Exp2.left.unaff.flip.flip.set'
'P16_Exp2.Right.unaff.flip.set'
'P18_Exp2.Right.unaff.flip.set'
'P20_Exp2.Right.Unaff.flip.set'
};

B={
'C:\\EEGdata\\P4_Exp2.Right.aff.flip.ref.set'
'C:\\EEGdata\\P7_Exp2.Left.aff.ref.set'
'C:\\EEGdata\\P12_Exp2.Left.aff.ref.set'
'C:\\EEGdata\\P15_Exp2_Right_aff.flip.ref.set'
'C:\\EEGdata\\P17_Exp2.Left.aff.ref.set'
'C:\\EEGdata\\P19_Exp2.Left.aff.ref.set'
'C:\\EEGdata\\P8_Exp2.Left.aff.ref.set'
'C:\\EEGdata\\P10_Exp2.right.aff.flip.ref.set'
'C:\\EEGdata\\P13_Exp2_Left.aff.ref.set'
'C:\\EEGdata\\P14_Exp2.right.aff.flip.ref.set'
'C:\\EEGdata\\P16_Exp2.Left.aff.ref.set'
'C:\\EEGdata\\P18_Exp2.Left.aff.ref.set'
'C:\\EEGdata\\P20_Exp2.Left.aff.ref.set'
'C:\\EEGdata\\P4_Exp2.Left.Unaff.flip.flip.ref.set'
'C:\\EEGdata\\P7_Exp2.Right.Unaff.flip.ref.set'
'C:\\EEGdata\\P12_Exp2.Right.Unaff.flip.ref.set'
'C:\\EEGdata\\P15_Exp2_Left_unaff.flip.flip.ref.set'
'C:\\EEGdata\\P17_Exp2.Right.unaff.flip.ref.set'
'C:\\EEGdata\\P19_Exp2.Right.Unaff.flip.ref.set'
'C:\\EEGdata\\P8_Exp2.Right.Unaff.flip.ref.set'
'C:\\EEGdata\\P10_Exp2.left.unaff.flip.flip.ref.set'
'C:\\EEGdata\\P13_Exp2_Right.unaff.flip.ref.set'
'C:\\EEGdata\\P14_Exp2.left.unaff.flip.flip.ref.set'
'C:\\EEGdata\\P16_Exp2.Right.unaff.flip.ref.set'
'C:\\EEGdata\\P18_Exp2.Right.unaff.flip.ref.set'
'C:\\EEGdata\\P20_Exp2.Right.Unaff.flip.ref.set'
};

subn=length(A);


for i=1:subn
    
[ALLEEG EEG CURRENTSET ALLCOM] = eeglab;
EEG = pop_loadset('filename',A(i),'filepath','C:\\EEGdata\\');
[ALLEEG, EEG, CURRENTSET] = eeg_store( ALLEEG, EEG, 0 );

EEG = eeg_checkset( EEG );
EEG = pop_reref( EEG, [35 81] );

[ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 1,'savenew',char(B(i)),'gui','off'); 

end