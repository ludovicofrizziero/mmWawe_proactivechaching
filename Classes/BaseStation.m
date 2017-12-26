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
        H %channel matrix
        H_params %stored channel matrix parameters
        f %carrier frequency
        BW %bandwidth
        AoA
        AoD
        %memory %for future use
    end
    
    properties
        ant_arr
        pos      
        BF %beam forming vector
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
        end
        
        function init(BS)
            %initialization
            BS_distance = norm(BS.sharedData.UE.pos - BS.pos);
            BS.PL = 10.^((32.4 + 21.6*log10(BS_distance) + 20*log10(BS.f/1e9))/10)';
            BS.find_AoA_AoD();
            [BS.H, BS.H_params] = compute_H_mobility(BS.f, BS.AoA, BS.AoD, BS.sharedData.UE.ant_arr, BS.ant_arr, BS.BW);
            
            BS.handin();
            if isequal(BS, BS.sharedData.servingBS)
                BS.BF = compute_BF_vector(BS.ant_arr, BS.AoA); %steer toward user
                BS.GAIN = (abs(conj(BS.sharedData.UE.BF) * (BS.sharedData.servingBS.H) * (BS.BF).').^2); %this station is transmitting
            else
                BS.BF = compute_BF_vector(BS.ant_arr, 2*pi*rand(1)); %random steering
                BS.GAIN = (abs(conj(BS.sharedData.UE.BF) * (BS.H) * (BS.BF).').^2); %this station is interfering
            end
                       
            BS.n = poissrnd(BS.mean_n);     
        end
        
        function update(BS, sim_time)                       
            if mod(sim_time, BS.th) < 1e-10
                %update all channel info
                BS.find_AoA_AoD();
                [BS.H, BS.H_params] = compute_H_mobility(BS.f, BS.AoA, BS.AoD, BS.sharedData.UE.ant_arr, BS.ant_arr, BS.BW);
            end
            
            if mod(sim_time, BS.tt) < 1e-10
                %update Beam Forming vector, keep channel params, update only ssf values
                BS_distance = norm(BS.sharedData.UE.pos - BS.pos);
                BS.PL = 10.^((32.4 + 21.6*log10(BS_distance) + 20*log10(BS.f/1e9))/10)';
                BS.find_AoA_AoD();
                BS.H = compute_H_ssf(BS.f, BS.AoA, BS.AoD, BS.sharedData.UE.ant_arr, BS.ant_arr, BS.H_params, BS.BW);
                
                BS.handin();
                if isequal(BS, BS.sharedData.servingBS)
                    BS.BF = compute_BF_vector(BS.ant_arr, BS.AoA); %steer toward user
                else
                    BS.BF = compute_BF_vector(BS.ant_arr, 2*pi*rand(1)); %random steering
                end
            else
                %don't update Beam Forming vector, keep channel params, update only ssf values
                BS.find_AoA_AoD();
                BS.H = compute_H_ssf(BS.f, BS.AoA, BS.AoD, BS.sharedData.UE.ant_arr, BS.ant_arr, BS.H_params, BS.BW);
            end
            
            if isequal(BS, BS.sharedData.servingBS)                
                BS.GAIN = (abs(conj(BS.sharedData.UE.BF) * (BS.sharedData.servingBS.H) * (BS.BF).').^2); %this station is transmitting
            else               
                BS.GAIN = (abs(conj(BS.sharedData.UE.BF) * (BS.H) * (BS.BF).').^2); %this station is interfering
            end
            
            if mod(sim_time, 1) < 1e10 
                %update every 1 sec the random n of users connected to this station
                BS.n = poissrnd(BS.mean_n); 
            end
        end
    end
    
    methods (Access = private)
        function find_AoA_AoD(BS)
            %AoA from the point of view of UE
            pos_x = -BS.sharedData.UE.pos(1) + BS.pos(1);
            pos_y = -BS.sharedData.UE.pos(2) + BS.pos(2);

            if (pos_x >= 0 && pos_y >=0)
                BS.AoD = atan(pos_y/pos_x);
            elseif (pos_x < 0 && pos_y >=0)
                BS.AoD = atan(pos_y/pos_x) + pi;
            elseif (pos_x < 0 && pos_y < 0)
                BS.AoD = atan(pos_y/pos_x) + pi;
            else
                BS.AoD = atan(pos_y/pos_x) + 2*pi;
            end
            
            BS.AoA = pi-BS.AoD;
        end
        
        function handin(BS)
            if BS.PL < BS.sharedData.servingBS.PL
                BS.sharedData.servingBS.handover(BS)
                BS.sharedData.servingBS = BS;
            end
            
        end
        
        function handover(BS, to_next_BS)
            %to_next_BS %reserved for future use
            BS.BF = compute_BF_vector(BS.ant_arr, 2*pi*rand(1)); %random steering, needed here to avoid errors
            BS.GAIN = (abs(conj(BS.sharedData.UE.BF) * (BS.H) * (BS.BF).').^2); %this station is interfering, needed here to avoid errors
        end
    end
    
end