%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Authors: Frizziero - Suman - Dell'Eva
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

classdef BaseStation < handle
    properties (Constant)
        sharedData = SharedBSData; %shared data among all BS instances
    end
    
    properties (SetAccess = private)         
        th
        tt
        mean_n
        n
        PL %path loss to the user
        C %channel state object
        f %carrier frequency
        BW %bandwidth
        AoA
        AoD
        ant_pos %all antenna positions
        max_memory
        memory
        TTT %Time To Transition
    end
    
    properties
        ID
        lambda %length of the carrier wave
        ant_arr
        pos      
        BF %beam forming object
        signal_power_at_ue
    end
    
    methods
        function BS = BaseStation(ID, antenna_array, f, BW, position, t_H, t_tracking, mean_nusers) %constructor
            BS.ID = ID;
            BS.ant_arr = antenna_array;
            BS.pos = position;
            BS.th = t_H;
            BS.tt = t_tracking;
            BS.TTT = 1; %for init
            BS.mean_n = mean_nusers;
            BS.f = f;
            BS.BW = BW;     
            BS.lambda = physconst('lightspeed') / f;
            BS.C = mmChannelModel(BW, f, 4, 2.8);            
            d = 0.15; % (m) spacing between antenna array elements 
            ind = 1;
            BS.ant_pos = zeros(BS.ant_arr, 3);
            for i = 1 : sqrt(BS.ant_arr)
                for j = 1 : sqrt(BS.ant_arr)
                    BS.ant_pos(ind, :) = [(i-1) * d , 0, (j-1) * d];
                    ind = ind + 1;
                end
            end          
            
            BS.BF = MVDR_Beamforming(0.01, BS.lambda, BS.ant_pos);
            BS.max_memory = int64(16e9 * 8); %[bits]
        end
        
        function init(BS)
            %initialization
            BS_distance = norm(BS.sharedData.UE.pos - BS.pos) / 1000;
            BS.PL = 10.^((32.4 + 20*log10(BS_distance) + 20*log10(BS.f/1e6) + 6*randn(1))/10);
            BS.find_AoD();            
            BS.C.update_channel_state(BS.sharedData.UE.AoA, BS.AoD, BS.sharedData.UE.ant_pos, BS.ant_pos, true);
            
            BS.signal_power_at_ue = 1/BS.PL; %just for init
            if isequal(BS, BS.sharedData.servingBS)                
                BS.BF.update_state(BS.AoD);     %this station is serving                      
            else        
                BS.BF.update_state([2*pi*rand(1), pi*rand(1)]); %this station is interfering
            end
            BS.compute_signal_power();
                       
            BS.n = poissrnd(BS.mean_n);
            M = 500e6 * 8; %mean
            V = 350e6 * 8; %variance
            used_mem = int64(sum(lognrnd(log(M^2/sqrt(V+M^2)), sqrt(log(V/M^2 + 1)), 1, BS.n)));            
            BS.memory = max(0, BS.max_memory - used_mem); %free memory
        end
        
        function update(BS, sim_time, dt)                            
            if mod(sim_time, BS.th) < 1e-10
                %update all channel info
                BS.find_AoD();
                BS.C.update_channel_state(BS.sharedData.UE.AoA, BS.AoD, BS.sharedData.UE.ant_pos, BS.ant_pos, true);
            end
            
            if mod(sim_time, BS.tt) < 1e-10
                %update Beam Forming vector, keep channel params, update only ssf values
                BS_distance = norm(BS.sharedData.UE.pos - BS.pos) / 1000; %Km
                BS.PL = 10^((32.4 + 20*log10(BS_distance) + 20*log10(BS.f/1e6) + 6*randn(1))/10);                                    
                BS.find_AoD();
                BS.C.update_channel_state(BS.sharedData.UE.AoA, BS.AoD, BS.sharedData.UE.ant_pos, BS.ant_pos, false);

                BS.handin(dt);
                if isequal(BS, BS.sharedData.servingBS) 
                    %% see DOC/meanconntime_bias.jpg for reference
                    a = pi/6; %30 deg                   
                    if (BS.AoD(1) >= a && BS.AoD(1) <= pi-a) || (BS.AoD(1) >= -(pi-a)  && BS.AoD(1) <= -a)
                        BS.BF.update_state(BS.AoD); %UE is in useful sweep range (+- 60 deg)
                    else
                        BS.BF.update_state([pi/2, BS.AoD(2)]); %UE is outside useful sweep range (+- 60 deg)
                    end
                else                    
                    BS.BF.update_state([2*pi*rand(1), pi*rand(1)]);
                end                
            else
                %don't update Beam Forming vector, keep channel params, update only ssf values
                BS.find_AoD();
                BS.C.update_channel_state(BS.sharedData.UE.AoA, BS.AoD, BS.sharedData.UE.ant_pos, BS.ant_pos, false);
            end
                                      
            BS.compute_signal_power(); %this station is interfering if it is not equal to sharedData.servingBS            
            
            if mod(sim_time, 1) < 1e10
                %update every 1 sec the random n of users connected to this station
                n_ = poissrnd(BS.mean_n); 
                if n_ > 1 
                    BS.n = n_;
                else
                    BS.n = 1;
                end
            end
        end
        
        function mem = get_mem_for_UE(BS)
            %at this point the station is offering a certain amount of space for the UE. 
            %It need not be the case the BS will actually receive the file
            %% see DOC/meanconntime_bias.jpg for reference
            mean_conn_time = ( abs(BS.pos(2) - BS.sharedData.UE.pos(2)) * sin(pi/3) / sin(pi/6) ) * 2 / BS.sharedData.UE.vel;
            %%
            mem = int64(BS.sharedData.UE.requested_rate * 10 * mean_conn_time * 1e9); %[bits]
            if mem > BS.memory
                mem = 0;
            end
        end
        
        function allocate_memory_for_ue(BS, mem)
            %at this point the station has been allocated with a chunck of the file requested by the UE
            BS.memory = mem; %free memory no more needed, reuse it as actual storage for the UE 
        end
        
        function file_chunk = download_file(BS, dt, rate)
            file_chunk = 0;
            if isequal(BS, BS.sharedData.servingBS)  
                if ~(BS.memory == 0)
                    file_chunk = int64(round(rate * dt));
                    BS.memory = BS.memory - file_chunk;
                    if BS.memory < 0
                        file_chunk = file_chunk + BS.memory;
                        BS.memory = 0;
                    end
                end
            end
        end
    end
    
    methods (Access = private)
        function find_AoD(BS)
            %for reference see DOC/BeamForming/angles.jpg
            
            %change point of view from the world origin to the BS's coord system
            new_ue_pos = BS.sharedData.UE.pos - BS.pos;     
            new_ue_pos = new_ue_pos / norm(new_ue_pos);
            
            theta = acos(dot(new_ue_pos(1:2), [1, 0])); %theta > 0 always, use only xy coords
            phi = acos(dot(new_ue_pos, [0,0,1])); %should be pi/2 <= phi < pi  because BS at least as tall as UE or more 
            tmp = cross([1, 0, 0], [new_ue_pos(1), new_ue_pos(2), 0]);
            s = sign(tmp(3));
            
            %find azimut angle (with respect to x axis, positive toward y axis)               
            if theta <= pi/2 && s > 0
                BS.AoD = [theta, phi];
            elseif theta < pi/2 && s < 0
                BS.AoD = [-theta, phi];
            elseif theta >= pi/2 && s > 0
                BS.AoD = [theta, phi];
            elseif theta > pi/2 && s < 0
                BS.AoD = [-theta, phi];
            end                      
        end               
        
        function handin(BS, dt)
%             if BS.PL < BS.sharedData.servingBS.PL
%                 BS.sharedData.servingBS.handover(BS)
%                 BS.sharedData.servingBS = BS;
%             end
            
            outage_thresh = 10^( -5 /10);
            thermal_noise = 10^((-174+7)/10) * (BS.BW / BS.n );
            SNR = BS.signal_power_at_ue / thermal_noise;
            if ~isequal(BS, BS.sharedData.servingBS) && (SNR > outage_thresh)
                minTTT = 0.3;  
                %% see DOC/meanconntime_bias.jpg for reference
                length = abs(BS.pos(2) - BS.sharedData.UE.pos(2)) * sin(pi/3) / sin(pi/6); %optimal  BS' cell **HALF** chord length relative to UE               
                delta = (BS.sharedData.UE.pos(1) - (BS.pos(1) - (length + 25)) )/ ((length + 25)*2 ); % -+25m tolerance
                bias = (BS.memory > 0) * delta * (delta >= 0 && delta <= 1);
                %%         
                P = BS.signal_power_at_ue * (1 + bias); %similar idea as presented in DOC/iswcs-symbiocity.pdf                 
                if  P > BS.sharedData.servingBS.signal_power_at_ue && BS.TTT >= minTTT % BS.PL < BS.sharedData.servingBS.PL && BS.TTT >= 0.5 % 
                    BS.sharedData.servingBS.handover(BS) %ask the old BS to hand over the UE
                    BS.sharedData.servingBS = BS;
                    BS.TTT = 0;
                elseif P > BS.sharedData.servingBS.signal_power_at_ue && BS.TTT < minTTT % BS.PL < BS.sharedData.servingBS.PL && BS.TTT < 0.5 % 
                    BS.TTT = BS.TTT + dt;
                else
                    BS.TTT = 0;
                end     
            end
        end
        
        function handover(BS, to_next_BS)
%             fprintf('handover from %d to %d;\n', BS.ID, to_next_BS.ID);
%             fprintf('UE xpos: %3.3f\t BS xpos: %3.3f\n', BS.sharedData.UE.pos(1), to_next_BS.pos(1));
            BS.BF.update_state([2*pi*rand(1), pi*rand(1)]);
            BS.compute_signal_power(); %this station is interfering, needed here to avoid errors
        end      
                
        function compute_signal_power(BS)
            P_tx = 0.5; %Transmitting power [Watt]
            BS.signal_power_at_ue = P_tx * abs((BS.sharedData.UE.BF)' * ((BS.C) * (BS.BF) + BaseStation.cnoise(BS.sharedData.UE.ant_arr, 0.01)) )^2 / BS.PL;
        end
    end
    
    methods (Static)
        function n = cnoise(N, var)
            %complex circular symmetric noise
            n = sqrt(var/2) * (randn(N,1) + 1i*randn(N,1));            
        end
    end
    
end