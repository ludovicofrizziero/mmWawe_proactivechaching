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
        ant_pos %all antenna positions in wawelength units
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
            BS.memory = round(rand(1));
        end
        
        function init(BS)
            %initialization
            BS_distance = norm(BS.sharedData.UE.pos - BS.pos) / 1000;
            BS.PL = 10.^((32.4 + 20*log10(BS_distance) + 20*log10(BS.f/1e6) + 6*randn(1))/10);
            BS.find_AoD();            
            BS.C.update_channel_state(BS.sharedData.UE.AoA, BS.AoD, BS.sharedData.UE.ant_pos, BS.ant_pos, true);
            
            BS.signal_power_at_ue = 1/BS.PL; %just for init
%             BS.handin();
            if isequal(BS, BS.sharedData.servingBS)                
                BS.BF.update_state(BS.AoD);     %this station is serving                      
            else        
                BS.BF.update_state([2*pi*rand(1), pi*rand(1)]); %this station is interfering
            end
            BS.compute_signal_power();
                       
            BS.n = poissrnd(BS.mean_n);     
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
                BS.n = poissrnd(BS.mean_n); 
            end
        end
        
        function download_file(BS, dt, rate)
            if isequal(BS, BS.sharedData.servingBS)
                %% download a portion (rate * dt) of the stored file
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
            
            P = BS.signal_power_at_ue; %similar idea as presented in DOC/iswcs-symbiocity.pdf but the bias is very different in nature
            if  P > BS.sharedData.servingBS.signal_power_at_ue && BS.TTT >= 0.2 % BS.PL < BS.sharedData.servingBS.PL && BS.TTT >= 0.5 % 
                BS.sharedData.servingBS.handover(BS) %ask the old BS to hand over the UE
                BS.sharedData.servingBS = BS;
                BS.TTT = 0;
            elseif P > BS.sharedData.servingBS.signal_power_at_ue && BS.TTT < 0.2 % BS.PL < BS.sharedData.servingBS.PL && BS.TTT < 0.5 % 
                BS.TTT = BS.TTT + dt;
            else
                BS.TTT = 0;
            end            
        end
        
        function handover(BS, to_next_BS)
            fprintf('handover from %d to %d;\n', BS.ID, to_next_BS.ID);
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