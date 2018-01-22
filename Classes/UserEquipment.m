%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Authors: Frizziero - Suman - Dell'Eva
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

classdef UserEquipment < handle
%     properties (Constant)
%         sharedData = SharedBSData
%     end
    
    properties (SetAccess = private)       
        tt
        lambda
        ant_pos
    end
    
    properties
        sharedData 
        pos  
        ant_arr
        BF %Beam Forming object
        vel
        AoD
        AoA
        requested_rate
        requested_file_size
        buffer
        max_buffer
        lost_data
    end
    
    methods   
        function UE = UserEquipment(antenna_array, f, position, vel, t_tracking) %constructor
            UE.ant_arr = antenna_array;
            UE.pos = position;
            UE.vel = vel;
            UE.tt = t_tracking;   
            UE.lambda = physconst('lightspeed') / f;            
            d = 0.01; %m
            UE.ant_pos = zeros(UE.ant_arr, 3);           
            ind = 1;
            for i = 1 : sqrt(UE.ant_arr)
                for j = 1 : sqrt(UE.ant_arr)
                    UE.ant_pos(ind, :) = [(i-1) * d , 0, (j-1) * d];
                    ind = ind + 1;
                end
            end
            
            UE.BF = MVDR_Beamforming(0.01, UE.lambda, UE.ant_pos);
            rates = [0.016, 0.024, 0.035, 0.045, 0.053, 0.068]; %1440p and 2160p (4K) common bitrates (compressed video) [Gbps]
            UE.requested_rate = rates(randi(size(rates), 1));
            M = 1e9 * 8;   %[bits] mean
            V = 350e6 * 8; %[bits] variance            
            UE.requested_file_size = int64(sum(lognrnd(log(M^2/sqrt(V+M^2)), sqrt(log(V/M^2 + 1))))); 
            UE.max_buffer = int64(1.5e9 * 8); %[bits]
            UE.buffer = 0;
        end        
        
        function init(UE)
            %init  
            UE.find_AoA()
            UE.BF.update_state(UE.AoA);
        end
        
        function update(UE, sim_time, dt)
            
            UE.pos(1) = dt * UE.vel + UE.pos(1); %+ UE.pos(1);
            
            if mod(sim_time, UE.tt) < 1e-10
                %update Beam Forming vector
                UE.find_AoA()
                UE.BF.update_state(UE.AoA);
            end
                       
            if UE.buffer > 0
                UE.buffer = UE.buffer - UE.requested_rate * dt;
                if UE.buffer < 0
                    %TODO signal bad QoS
                    UE.lost_data = UE.lost_data - UE.buffer;
                    
                    UE.buffer = 0;
                end
            end
        end
        
        function receive_file_chunk(UE, file_chunk)
            UE.buffer = UE.buffer + file_chunk;
            if UE.buffer > UE.max_buffer
                %TODO signal bad QoS
                UE.lost_data = UE.lost_data + (UE.buffer - UE.max_buffer);
                
                UE.buffer = UE.max_buffer; 
            end
        end
    end
    
    methods (Access = private)
        function find_AoA(UE)
            %for reference see DOC/BeamForming/angles.jpg
            
            %change point of view from the world origin to the BS's system
            new_ue_pos = UE.pos - UE.sharedData.servingBS.pos;    
            new_ue_pos = new_ue_pos / norm(new_ue_pos);
            
            theta = acos(dot(-new_ue_pos(1:2), [1, 0])); %only xy coords 
            phi = acos(dot(-new_ue_pos, [0,0,1])); %should be 0 < phi <= pi/2 because BS is as tall as UE or more
            tmp = cross([1, 0, 0], [new_ue_pos(1), new_ue_pos(2), 0]);
            s = sign(tmp(3));
            
            %find azimut angle (with respect to y axis, positive toward x axis)            
            if theta <= pi/2 && s > 0 
                UE.AoA = [theta, phi];
            elseif theta < pi/2 && s < 0
                UE.AoA = [-theta, phi];
            elseif theta >= pi/2 && s > 0
                UE.AoA = [theta, phi];
            elseif theta > pi/2 && s < 0
                UE.AoA = [-theta, phi];
            end           
        end
    end
end