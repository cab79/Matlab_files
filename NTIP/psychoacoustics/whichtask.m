function whichtask (sel,logicswp)
global h
mh = guihandles(gcf);
switch sel
    case 'yes/no'
        h.tasktype = 0;
        set(mh.textcatchtrial,'Visible', 'on');
        set(mh.catchprop,'Visible', 'on');
        set(mh.nafc,'Visible', 'off');
        set(mh.textnafc,'Visible','off');
        set(mh.nafc,'String', 0);
        set(mh.feedback,'Value', 0);
        if logicswp
            set(mh.p_target,'String', '0.631');
            set (mh.swp,'Value',1);
        else
            set(mh.p_target,'String', '0.5');
            set (mh.swp,'Value',0);
        end
    case 'nAFC'
        h.tasktype = 1;
        nafc= str2double(get(mh.nafc,'String'));
        if nafc
            set(mh.nafc,'String',nafc);
        else
            set(mh.nafc,'String',2);
        end
        set(mh.ptargetupdate,'Visible', 'on');
        gamma = 1/str2double( get(mh.nafc,'String'));
        set(mh.nafc,'Visible', 'on');
        set(mh.textnafc,'Visible','on');
        set(mh.textcatchtrial,'Visible', 'of');
        set(mh.catchprop,'Visible', 'of');

        if logicswp
            swp = (2*gamma+1+sqrt(1+8*gamma))/(3+sqrt(1+8*gamma));
            set(mh.p_target,'String', swp);
            set (mh.swp,'Value',1);
        else
            set(mh.p_target,'String', (1+gamma)/2);
            set (mh.swp,'Value',0);
        end
end
