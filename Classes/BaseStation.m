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
        AoA_compatibility %for compatibility with old code
        AoD_compatibility %for compatibility with old code
        ant_pos %all antenna positions in wawelength units
        %memory %for future use
    end
    
    properties
        lambda %length of the carrier wawe
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
            BS.lambda = physconst('lightspeed') / f;
            d = 0.15; % (m) spacing between antenna array elements 
            ind = 1;
            BS.ant_pos = zeros(BS.ant_arr, 3);
            for i = 1 : sqrt(BS.ant_arr)
                for j = 1 : sqrt(BS.ant_arr)
                    BS.ant_pos(ind, :) = [(i-1) * d , 0, (j-1) * d];
                    ind = ind + 1;
                end
            end          
            
        end
        
        function init(BS)
            %initialization
            BS_distance = norm(BS.sharedData.UE.pos - BS.pos);
            BS.PL = 10.^((32.4 + 21.6*log10(BS_distance) + 20*log10(BS.f/1e9))/10)';
            BS.find_AoA_AoD();
            [BS.H, BS.H_params] = compute_H_mobility(BS.f, BS.sharedData.UE.AoA, BS.AoD, BS.sharedData.UE.ant_pos, BS.ant_pos, BS.BW);
            
            BS.handin();
            if isequal(BS, BS.sharedData.servingBS)
                %BS.BF = compute_BF_vector(BS.ant_arr, BS.AoA); %steer toward user
                BS.BF = compute_BF_vector2(BS.AoD, BS.lambda, BS.ant_pos);               
                BS.GAIN = (abs(conj(BS.sharedData.UE.BF).' * (BS.sharedData.servingBS.H) * (BS.BF)).^2); %this station is transmitting
            else
                %BS.BF = compute_BF_vector(BS.ant_arr, 2*pi*rand(1)); %random steering
                BS.BF = compute_BF_vector2([2*pi*rand(1), 2*pi*rand(1)], BS.lambda, BS.ant_pos);
                BS.GAIN = (abs(conj(BS.sharedData.UE.BF).' * (BS.H) * (BS.BF)).^2); %this station is interfering
            end
                       
            BS.n = poissrnd(BS.mean_n);     
        end
        
        function update(BS, sim_time)                       
            if mod(sim_time, BS.th) < 1e-10
                %update all channel info
                BS.find_AoA_AoD();
                [BS.H, BS.H_params] = compute_H_mobility(BS.f, BS.sharedData.UE.AoA, BS.AoD, BS.sharedData.UE.ant_pos, BS.ant_pos, BS.BW);
            end
            
            if mod(sim_time, BS.tt) < 1e-10
                %update Beam Forming vector, keep channel params, update only ssf values
                BS_distance = norm(BS.sharedData.UE.pos - BS.pos) / 1000; %Km
                BS.PL = 10.^((32.4 + 21.6*log10(BS_distance) + 20*log10(BS.f/1e9))/10);
                BS.find_AoA_AoD();
                BS.H = compute_H_ssf(BS.f, BS.sharedData.UE.AoA, BS.AoD, BS.sharedData.UE.ant_pos, BS.ant_pos, BS.H_params, BS.BW);
                
                BS.handin();
                if isequal(BS, BS.sharedData.servingBS)
                    %BS.BF = compute_BF_vector(BS.ant_arr, BS.AoA); %steer toward user
                    BS.BF = compute_BF_vector2(BS.AoD, BS.lambda, BS.ant_pos);
                else
                    %BS.BF = compute_BF_vector(BS.ant_arr, 2*pi*rand(1)); %random steering
                    BS.BF = compute_BF_vector2([2*pi*rand(1), 2*pi*rand(1)], BS.lambda, BS.ant_pos);
                end
            else
                %don't update Beam Forming vector, keep channel params, update only ssf values
                BS.find_AoA_AoD();
                BS.H = compute_H_ssf(BS.f, BS.sharedData.UE.AoA, BS.AoD, BS.sharedData.UE.ant_pos, BS.ant_pos, BS.H_params, BS.BW);
            end
                                      
            BS.GAIN = (abs((BS.sharedData.UE.BF)' * (BS.H) * (BS.BF)).^2); %this station is interfering if it is not equal to sharedData.servingBS
            
%             %for debug
%             if isequal(BS, BS.sharedData.servingBS) && BS.ant_arr > 16 
%                 disp('----')
%                 disp(abs(BS.sharedData.UE.BF' * BS.H * BS.BF)^2)
%                 disp('----')
%             elseif BS.ant_arr > 16
%                 disp(abs(BS.sharedData.UE.BF' * BS.H * BS.BF)^2)
%             end
            
            if mod(sim_time, 1) < 1e10
                %update every 1 sec the random n of users connected to this station
                BS.n = poissrnd(BS.mean_n); 
            end
        end
    end
    
    methods (Access = private)
        function find_AoA_AoD(BS)
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
            %to_next_BS %reserved for future use
            %BS.BF = compute_BF_vector(BS.ant_arr, 2*pi*rand(1)); %random steering, needed here to avoid errors
            BS.BF = compute_BF_vector2([2*pi*rand(1), 2*pi*rand(1)], BS.lambda, BS.ant_pos);
            BS.GAIN = (abs(conj(BS.sharedData.UE.BF).' * (BS.H) * (BS.BF)).^2); %this station is interfering, needed here to avoid errors
        end
    end
    
end