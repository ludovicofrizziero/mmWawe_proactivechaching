function distance_handover(allBS, shared_data)
    UE = shared_data.UE;
    if shared_data.servingBS.memory <= 0
        d = zeros(length(allBS), 1);
        for i = 1:length(allBS)
            tmp = abs(allBS{i}.pos(1) - UE.pos(1));
            if allBS{i}.memory > 0 && UE.pos(1) >= allBS{i}.support_start
                d(i) = tmp;
            else
                d(i) = 1e10;
            end
        end
        [~, i] = min(d);

        shared_data.servingBS.handover(allBS{i});
    end
end