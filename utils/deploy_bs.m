function [n_BS_top, n_BS_bottom, pos_top, pos_bottom] = deploy_bs(BS_per_km, road_length, R)
    n_BS_top = poissrnd(BS_per_km/1000 * (road_length),1,1); % number of BSs
    n_BS_bottom = poissrnd(BS_per_km/1000 * (road_length),1,1); % number of BSs   
%     n_BS_top = BS_per_km;
%     n_BS_bottom = BS_per_km;
    BS_distance_avg_top = 1000/n_BS_top;
    BS_distance_avg_bottom = 1000/n_BS_bottom;
%     delta_top = ceil(BS_distance_avg_top/8);
%     delta_bottom = ceil(BS_distance_avg_bottom/8); 
    delta_top = 500/BS_per_km;
    delta_bottom = 500/BS_per_km; 


    pos_top = zeros(n_BS_top, 3);
    for i = 1:n_BS_top
        pos_top(i, :) = [(i-1)*BS_distance_avg_top+unifrnd(-delta_top , delta_top) , 2 * R, 8];            
    end

    pos_bottom = zeros(n_BS_bottom, 3);
    for i = 1:n_BS_bottom
        pos_bottom(i, :) = [500/BS_per_km * 0 + (i-1)*BS_distance_avg_bottom+unifrnd(-delta_bottom , delta_bottom) , 0, 8];
    end    

%     disp(n_BS_bottom + n_BS_top);
end