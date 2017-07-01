% LOR2SPM Batch

clear all

subjects = {'H1', 'H2', 'H3', 'H4', 'H5', 'H6', 'H7', 'H8', 'H9', 'H10', 'H11', 'H12', 'H13' 'H14', 'H15', 'H16', 'F1', 'F2', 'F3', 'F4', 'F5', 'F6', 'F7', 'F8', 'F9', 'F10', 'F11', 'F12', 'F13', 'F14', 'F15', 'F16', 'OA1', 'OA2', 'OA4', 'OA5', 'OA6', 'OA7', 'OA8', 'OA9', 'OA10', 'OA11', 'OA12','OA13','OA14', 'OA15', 'OA16','OA17'};
   
Nsub = length(subjects);



 fnames={'1_e_av-B2T.lor';
    '2_e_av-B2T.lor';
    '3_e_av-B2T.lor';
    '4_e_av-B2T.lor';
    '5_e_av-B2T.lor';
    '6_e_av-B2T.lor';
    '1_l_av-B2T.lor';
    '2_l_av-B2T.lor';
    '3_l_av-B2T.lor';
    '4_l_av-B2T.lor';
    '5_l_av-B2T.lor';
    '6_l_av-B2T.lor';
    '1_m_av-B2T.lor';
    '2_m_av-B2T.lor';
    '3_m_av-B2T.lor';
    '4_m_av-B2T.lor';
    '5_m_av-B2T.lor';
    '6_m_av-B2T.lor';
    '1_p2_av-B2T.lor';
    '2_p2_av-B2T.lor';
    '3_p2_av-B2T.lor';
    '4_p2_av-B2T.lor';
    '5_p2_av-B2T.lor';
    '6_p2_av-B2T.lor';
    '1_n2_av-B2T.lor';
    '2_n2_av-B2T.lor';
    '3_n2_av-B2T.lor';
    '4_n2_av-B2T.lor';
    '5_n2_av-B2T.lor';
    '6_n2_av-B2T.lor';
    };




for n = 1:Nsub
    subject = subjects(n);
    subject = char(subject);   

    for x = 1:length(fnames)
    
    fname=char(fnames(x));
    fname2= [subject '_' fname];
    l=load (fname2);
    %l=1000*l;
    l=log(l);
    eval (['save ' subject '_log10_' fname ' l /ASCII'])
    
    end
   
end
    




