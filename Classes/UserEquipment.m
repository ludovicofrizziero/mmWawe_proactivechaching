classdef UserEquipment < handle
%     properties (Constant)
%         sharedData = SharedBSData
%     end
    
    properties (SetAccess = private)
        AoD
        tt
    end
    
    properties
        sharedData 
        pos  
        ant_arr
        BF %Beam Forming vector
        vel
    end
    
    methods   
        function UE = UserEquipment(antenna_array, position, vel, t_tracking) %constructor
            UE.ant_arr = antenna_array;
            UE.pos = position;
            UE.vel = vel;
            UE.tt = t_tracking;                        
        end        
        
        function init(UE)
            %init  
            UE.find_AoD()
            UE.BF = compute_BF_vector(UE.ant_arr, UE.AoD);
        end
        
        function update(UE, sim_time)
            
            UE.pos(1) = sim_time * UE.vel; %+ UE.pos(1);
            
            if mod(sim_time, UE.tt) < 1e-10
                %update Beam Forming vector
                UE.find_AoD()
                UE.BF = compute_BF_vector(UE.ant_arr, UE.AoD);                                       
            end
            
        end
    end
    
    methods (Access = private)
        function find_AoD(UE)
            %AoD from the point of view of UE
            pos_x = -UE.pos(1) + UE.sharedData.servingBS.pos(1);
            pos_y = -UE.pos(2) + UE.sharedData.servingBS.pos(2);

            if (pos_x >= 0 && pos_y >=0)
                UE.AoD = atan(pos_y/pos_x);
            elseif (pos_x < 0 && pos_y >=0)
                UE.AoD = atan(pos_y/pos_x) + pi;
            elseif (pos_x < 0 && pos_y < 0)
                UE.AoD = atan(pos_y/pos_x) + pi;
            else
                UE.AoD = atan(pos_y/pos_x) + 2*pi;
            end
        end
    end
end