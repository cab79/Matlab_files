function ForwardMaskingPar_MLP(mh)
        set(mh.beta,'String', '1');
        set(mh.nblocks,'String', '5');
        set(mh.ntrials,'String', '30');
        whichtask('nAFC',1);
        set(mh.nafcbutton,'value', 1);
        set(mh.nafc, 'String', '3');
        set(mh.hypnr,'String', '100');
        set(mh.firstmidpoint,'String', '-50');
        set(mh.lastmidpoint,'String', '-5');
        set(mh.startinglevel,'String', '-10');
        set(mh.logscale, 'value', 0);