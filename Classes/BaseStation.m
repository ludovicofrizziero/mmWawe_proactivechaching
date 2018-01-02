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
        AoA_compatibility %for compatibility with old code
        AoD_compatibility %for compatibility with old code
        ant_pos %all antenna positions in wawelength units
        %memory %for future use
    end
    
    properties
        lambda %length of the carrier wawe
        ant_arr
        pos      
        BF %beam forming object
        GAIN
    end
    
    methods
        function BS = BaseStation(antenna_array, f, BW, position, t_H, t_tracking, mean_nusers) %constructor
            BS.ant_arr = antenna_array;
            BS.pos = position;
            BS.th = t_H;
            BS.tt = t_tracking;
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
        end
        
        function init(BS)
            %initialization
            BS_distance = norm(BS.sharedData.UE.pos - BS.pos) / 1000;
            %BS.PL = 10.^((32.4 + 20*log10(BS_distance) + 20*log10(BS.f/1e9))/10);
            BS.PL = 10^(22.7 + 36.7 * log10(BS_distance * 1000) + 26 * log10(BS.f / 1e9) / 10);
            BS.find_AoD();            
            BS.C.update_channel_state(BS.sharedData.UE.AoA, BS.AoD, BS.sharedData.UE.ant_pos, BS.ant_pos, true);
            
            BS.handin();
            if isequal(BS, BS.sharedData.servingBS)                
                BS.BF.update_state(BS.AoD);               
                BS.GAIN = (abs((BS.sharedData.UE.BF)' * (BS.sharedData.servingBS.C) * (BS.BF)).^2); %this station is transmitting
            else        
                BS.BF.update_state([2*pi*rand(1), 2*pi*rand(1)]);
                BS.GAIN = (abs((BS.sharedData.UE.BF)' * (BS.C) * (BS.BF)).^2); %this station is interfering
            end
                       
            BS.n = poissrnd(BS.mean_n);     
        end
        
        function update(BS, sim_time)                       
            if mod(sim_time, BS.th) < 1e-10
                %update all channel info
                BS.find_AoD();
                BS.C.update_channel_state(BS.sharedData.UE.AoA, BS.AoD, BS.sharedData.UE.ant_pos, BS.ant_pos, true);
            end
            
            if mod(sim_time, BS.tt) < 1e-10
                %update Beam Forming vector, keep channel params, update only ssf values
                BS_distance = norm(BS.sharedData.UE.pos - BS.pos) / 1000; %Km
                %BS.PL = 10.^((32.4 + 20*log10(BS_distance) + 20*log10(BS.f/1e9))/10);
                BS.PL = 10^(22.7 + 36.7 * log10(BS_distance * 1000) + 26 * log10(BS.f / 1e9) / 10);
                BS.find_AoD();
                BS.C.update_channel_state(BS.sharedData.UE.AoA, BS.AoD, BS.sharedData.UE.ant_pos, BS.ant_pos, false);

                BS.handin();
                if isequal(BS, BS.sharedData.servingBS)                    
                    BS.BF.update_state(BS.AoD);
                else                    
                    BS.BF.update_state([2*pi*rand(1), 2*pi*rand(1)]);
                end
            else
                %don't update Beam Forming vector, keep channel params, update only ssf values
                BS.find_AoD();
                BS.C.update_channel_state(BS.sharedData.UE.AoA, BS.AoD, BS.sharedData.UE.ant_pos, BS.ant_pos, false);
            end
                                      
            BS.GAIN = (abs((BS.sharedData.UE.BF)' * (BS.C) * (BS.BF)).^2); %this station is interfering if it is not equal to sharedData.servingBS            
            
            if mod(sim_time, 1) < 1e10
                %update every 1 sec the random n of users connected to this station
                BS.n = poissrnd(BS.mean_n); 
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
        
        function handin(BS)
            if BS.PL < BS.sharedData.servingBS.PL
                BS.sharedData.servingBS.handover(BS)
                BS.sharedData.servingBS = BS;
            end
            
        end
        
        function handover(BS, to_next_BS)
            BS.BF.update_state([2*pi*rand(1), 2*pi*rand(1)]);
            BS.GAIN = (abs((BS.sharedData.UE.BF)' * (BS.C) * (BS.BF)).^2); %this station is interfering, needed here to avoid errors
        end
    end
    
end