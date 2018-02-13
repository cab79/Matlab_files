function open_serial(h,opt)

switch opt
    case 'spt1'
        if exist('spt1','var')
            try
                delete(spt1);
            end
        end
        global spt1
        PortNum=h.Settings.spt1_port; % PortNum variable is a string use eg: 'LPT1'
        spt1 = serial(PortNum); 
        fopen(spt1);
    case 'spt'
        global spt
        spt = serial(h.Settings.spt_port); 
        set(spt,'BaudRate',9600);
        set(spt,'DataBits',8);
        set(spt,'FlowControl','none');
        set(spt,'Parity','none');
        set(spt,'StopBits',1);
        set(spt,'Terminator','');
        fopen(spt);
end