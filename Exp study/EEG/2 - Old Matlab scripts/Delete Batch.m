% Delete Batch

clear all

% LOR2SPM Batch

clear all

subjects = {'H1', 'H2', 'H3', 'H4', 'H5', 'H6', 'H7', 'H8', 'H9', 'H10', 'H11', 'H12', 'H13' 'H14', 'H15', 'H16', 'F1', 'F2', 'F3', 'F4', 'F5', 'F6', 'F7', 'F8', 'F9', 'F10', 'F11', 'F12', 'F13', 'F14', 'F15', 'F16', 'OA1', 'OA2', 'OA4', 'OA5', 'OA6', 'OA7', 'OA8', 'OA9', 'OA10', 'OA11', 'OA12'};
  
Nsub = length(subjects);


 fnames={'1_e_av-B2T_F0001.img';
    '2_e_av-B2T_F0001.img';
    '3_e_av-B2T_F0001.img';
    '4_e_av-B2T_F0001.img';
    '5_e_av-B2T_F0001.img';
    '6_e_av-B2T_F0001.img';
    '1_l_av-B2T_F0001.img';
    '2_l_av-B2T_F0001.img';
    '3_l_av-B2T_F0001.img';
    '4_l_av-B2T_F0001.img';
    '5_l_av-B2T_F0001.img';
    '6_l_av-B2T_F0001.img';
    '1_m_av-B2T_F0001.img';
    '2_m_av-B2T_F0001.img';
    '3_m_av-B2T_F0001.img';
    '4_m_av-B2T_F0001.img';
    '5_m_av-B2T_F0001.img';
    '6_m_av-B2T_F0001.img';
    '1_p2_av-B2T_F0001.img';
    '2_p2_av-B2T_F0001.img';
    '3_p2_av-B2T_F0001.img';
    '4_p2_av-B2T_F0001.img';
    '5_p2_av-B2T_F0001.img';
    '6_p2_av-B2T_F0001.img';
    '1_n2_av-B2T_F0001.img';
    '2_n2_av-B2T_F0001.img';
    '3_n2_av-B2T_F0001.img';
    '4_n2_av-B2T_F0001.img';
    '5_n2_av-B2T_F0001.img';
    '6_n2_av-B2T_F0001.img';
    '1_vn_av-B2T_F0001.img';
    '2_vn_av-B2T_F0001.img';
    '3_vn_av-B2T_F0001.img';
    '4_vn_av-B2T_F0001.img';
    '5_vn_av-B2T_F0001.img';
    '6_vn_av-B2T_F0001.img';
    '1_vp_av-B2T_F0001.img';
    '2_vp_av-B2T_F0001.img';
    '3_vp_av-B2T_F0001.img';
    '4_vp_av-B2T_F0001.img';
    '5_vp_av-B2T_F0001.img';
    '6_vp_av-B2T_F0001.img';
    '1_e_av-B2T_F0001.hdr';
    '2_e_av-B2T_F0001.hdr';
    '3_e_av-B2T_F0001.hdr';
    '4_e_av-B2T_F0001.hdr';
    '5_e_av-B2T_F0001.hdr';
    '6_e_av-B2T_F0001.hdr';
    '1_l_av-B2T_F0001.hdr';
    '2_l_av-B2T_F0001.hdr';
    '3_l_av-B2T_F0001.hdr';
    '4_l_av-B2T_F0001.hdr';
    '5_l_av-B2T_F0001.hdr';
    '6_l_av-B2T_F0001.hdr';
    '1_m_av-B2T_F0001.hdr';
    '2_m_av-B2T_F0001.hdr';
    '3_m_av-B2T_F0001.hdr';
    '4_m_av-B2T_F0001.hdr';
    '5_m_av-B2T_F0001.hdr';
    '6_m_av-B2T_F0001.hdr';
    '1_p2_av-B2T_F0001.hdr';
    '2_p2_av-B2T_F0001.hdr';
    '3_p2_av-B2T_F0001.hdr';
    '4_p2_av-B2T_F0001.hdr';
    '5_p2_av-B2T_F0001.hdr';
    '6_p2_av-B2T_F0001.hdr';
    '1_n2_av-B2T_F0001.hdr';
    '2_n2_av-B2T_F0001.hdr';
    '3_n2_av-B2T_F0001.hdr';
    '4_n2_av-B2T_F0001.hdr';
    '5_n2_av-B2T_F0001.hdr';
    '6_n2_av-B2T_F0001.hdr';
    '1_vn_av-B2T_F0001.hdr';
    '2_vn_av-B2T_F0001.hdr';
    '3_vn_av-B2T_F0001.hdr';
    '4_vn_av-B2T_F0001.hdr';
    '5_vn_av-B2T_F0001.hdr';
    '6_vn_av-B2T_F0001.hdr';
    '1_vp_av-B2T_F0001.hdr';
    '2_vp_av-B2T_F0001.hdr';
    '3_vp_av-B2T_F0001.hdr';
    '4_vp_av-B2T_F0001.hdr';
    '5_vp_av-B2T_F0001.hdr';
    '6_vp_av-B2T_F0001.hdr';
    };




for n = 1:Nsub
    subject = subjects(n);
    subject = char(subject);   

    for x = 1:length(fnames)
    
    fname=char(fnames(x));
    fname2= [subject '_' fname];
    delete(fname2);
    
    end
   
end

