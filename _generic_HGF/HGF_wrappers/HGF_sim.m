function [D,S] = HGF_sim(D,S)

for d = 1:length(D) 

    %sim=[];
    for ns = 1:S.numsimrep
        S.HGF.selectrep=1;
        D_out=HGF_run(D(d),S,1);
        D(d).HGF(ns).sim = D_out.HGF(1).sim;
    end

    for ns = 1:S.numsimrep
        if length(D(d).HGF(ns).sim.y)==length(D(d).HGF(1).u)
            D(d).Output(ns).presstrial = 1:length(D(d).HGF(ns).sim.y);
            try 
                buttonopt = D(d).Output(1).Settings.buttonopt;
            catch
                buttonopt={};
            end
            if length(unique(D(d).HGF(ns).sim.y+1))==2 % should be binary, otherwise may be RT
                if ~isempty(buttonopt)
                    D(d).Output(ns).pressbutton = buttonopt(D(d).HGF(ns).sim.y+1);
                else
                    D(d).Output(ns).pressbutton = D(d).HGF(ns).sim.y;
                end
            end
        else
            error('simulated data has the wrong number of trials')
        end
    end
end