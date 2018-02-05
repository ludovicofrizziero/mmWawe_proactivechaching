function power_handover(allBS, shared_data)
    outage_thresh = 10^( -100 /10);
    N = max(size((allBS)));
    if shared_data.servingBS.memory <= 0       
        P = zeros(N, 1);
        for i = 1:N
            if allBS{i}.memory > 0
                thermal_noise = 10^((-174+7)/10) * (allBS{i}.BW / allBS{i}.n );
                SNR = allBS{i}.signal_power_at_ue / thermal_noise;
                if SNR > outage_thresh
                    bias = (1/abs(allBS{i}.pos(1) - shared_data.UE.pos(1)));
                    P(i) = (allBS{i}.memory > 0) * allBS{i}.signal_power_at_ue * bias;
                end
            end
        end
        
        [~, i] = max(P);        
        shared_data.servingBS.handover(allBS{i});
    end      

end