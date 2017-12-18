%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% V2I in vehicular network operating at millimeter wave.
% A multiu-lane highway scenario is considered.
%
% Marco Giordani
% July 2017
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all;
%close all;
clc;
rng(4)

print_scenario = false;
addpath(genpath('fnc'))

%% parameters for channel
f = 28e9; %carrier frequency [Hz]
c = 3e8; % light speed [m/s]
C_L = 10^(-7.2);
C_N = 10^(-6.14);
alpha_L = 2; % LOS PL exponent
alpha_N = 2.92; % NLOS PL exponent
P_tx = 27; % transmitting power [dBm]
P_tx_lin = 0.5;
BW = 1e9; % Bandwidth [Hz]
thermal_noise = 10^((-174+7)/10)*BW ; % thermal noise in linear scale (OBS: -30 to convert in dB)
sigma = thermal_noise/P_tx_lin; % thermal noise, normalized at the tranmitting power


n_tx = [4 16 64];
n_rx = [4  4 16];

T_sim = 20; % simulation duration [s]
t_offset = 0.1; %  simulation step [s]
T_tracking = 50; % tracking periodicity [s]
v = 50/3.6; % speed [m/s]
t_H = 0.3; % udpate of channel instances

road_length = 1000; % positive length of the road [m]
road_start = -1000; % negative length of the road [m]


%% parameters for simulation
n_rep_PL = 10;
theta_out = -5; %SINR outage threshold [dB]
theta_out_lin = 10.^(theta_out./10); %SINR outage threshold

%% paramers for road
d_R = 2.5; 
W_L = 3.5; % truck length
N_0 = 2; % # of obstacle lanes
R = d_R + N_0*W_L; % road width

UE = [road_length/2 R+ W_L/2]; % position of the UE

tic;
for antenna_idx = 1:length(n_tx) % consider single values of BS density
    for vector_lambda_bs = 30/1000;  
        for n_iter = 1:n_rep_PL
            
            fprintf('iteration: %d/%d -> %d/%d\n', antenna_idx, length(n_tx),n_iter, n_rep_PL);
            
            BS_all = [];
            %% deployment TOP part of road
            obs_lanes_top = cell(1,N_0);
   
            n_BS_top = poissrnd(vector_lambda_bs * (road_length),1,1); % number of BSs
            BS_top = [unifrnd(0, road_length, n_BS_top,1) , 2*R * ones(n_BS_top,1)]; %position of BSs
            
            n_BS_bottom = poissrnd(vector_lambda_bs * (road_length),1,1); % number of BSs
            BS_bottom = [unifrnd(0, road_length, n_BS_bottom,1) , 0 * ones(n_BS_bottom,1)]; %position of BSs
            
            BS_all = [BS_top; BS_bottom];
         
            n_BS_all = size(BS_all,1);
            BS_distance = abs( sqrt( (BS_all(:,1) - UE(1)).^2 + (BS_all(:,2) - UE(2)).^2 ) ); % new distance from BSs to current position of the UE
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %% PRINT SCENARIO (if enabled)
            if print_scenario
                hold on
                scatter(BS_all(:,1),BS_all(:,2),'*r'); hold on;
                scatter(UE(1),UE(2),'dm')
            end
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            %% Find PL and PL states, according to the distance
            BS_PL = zeros(1,n_BS_all);
            BS_PL = 10.^((32.4 + 21.6*log10(BS_distance) + 20*log10(f/1e9))/10)';
         
            %% MAKE ATTACHMENT DECISION (max. imparing PL, accoridng to the PL defined in Eq. (2))
            x_BS = BS_all(1,:); % x-coordinates of all BSs
            [serving_PL, serving_BS_idx] = min(BS_PL); % accoring to Assumption 3.4 and the definition of impariment PL in Eq. (2)
            distance_serving_BS = abs(BS_distance(serving_BS_idx)); % distance from UE in (0,0) to current serving BS
            
            %% AoA, AoD, spatial signatures, H computation for ues with respect to all BSs in the system
            AoA = zeros(n_BS_all,1);
            AoD = zeros(n_BS_all,1);
            H = cell(n_BS_all,1);
            H_params = cell(n_BS_all,1);
            for bs_idx = 1:n_BS_all
                pos_bs = BS_all(bs_idx,:);
                [tmp_AoA tmp_AoD] = find_angles(UE,pos_bs);
                AoA(bs_idx) = tmp_AoA; % AoA of 'ue' with respect to 'bs'
                AoD(bs_idx) = tmp_AoD;
                [H{bs_idx} H_params{bs_idx}] = compute_H_mobility(f,AoA(bs_idx),AoD(bs_idx),n_rx(antenna_idx),n_tx(antenna_idx)); % channel 
            end
            
            %% All BSs (but the serving one) steer in random directions, while the UE is aligned with the serving BS
            BF_vector_tx = cell(n_BS_all,1);
            rnd_steering_angle = 2*pi*rand(n_BS_all,1); % random steering angles
            for bs_idx = 1:n_BS_all
                BF_vector_tx{bs_idx} = compute_BF_vector(n_tx(antenna_idx),rnd_steering_angle(bs_idx));
            end
            BF_vector_tx{serving_BS_idx} = compute_BF_vector(n_tx(antenna_idx),AoA(serving_BS_idx));
            
            BF_vector_rx = compute_BF_vector(n_rx(antenna_idx),AoD(serving_BS_idx));
            Gain_num = (abs(conj(BF_vector_rx)*H{serving_BS_idx}*BF_vector_tx{serving_BS_idx}.').^2); 
            
            Gain_interference = zeros(n_BS_all,1);
            for bs_idx = 1:n_BS_all
                Gain_interference(bs_idx) = (abs(conj(BF_vector_rx)*H{bs_idx}*BF_vector_tx{bs_idx}.').^2);
            end
            
            SINR_interference = Gain_interference( setdiff(1:n_BS_all , serving_BS_idx),: ) ./ ...
                (BS_PL(setdiff(1:n_BS_all , serving_BS_idx))').*0.5;
            
            SINR_num =    Gain_num ./ BS_PL(serving_BS_idx)  ; % numerator of SINR (depends on the beamwidth)
            
            SINR_den = sum(SINR_interference) + thermal_noise;
            original_avg_SINR = SINR_num ./ SINR_den;
            original_rate = BW * log2(1+original_avg_SINR);
            
            %% consider time-varying simulation
            ue_test_new_position = [UE(1) UE(2)];
            index_internal = 1; %internal cycle index
            for t = t_offset:t_offset:T_sim
                %t
                space_offset = t_offset * v;
                ue_test_new_position = [ue_test_new_position(1)+space_offset ue_test_new_position(2)]; %new position of the user, after it has moved
                
                if mod(t, t_H) < 1e-10 %% change PL status for all BSs
                    BS_distance = abs( sqrt( (BS_all(:,1) - ue_test_new_position(1)).^2 + (BS_all(:,2) - ue_test_new_position(2)).^2 ) ); % new distance from BSs to current position of the UE
                    BS_PL = zeros(1,n_BS_all);
                    BS_PL = 10.^((32.4 + 21.6*log10(BS_distance) + 20*log10(f/1e9))/10)';
                    %% AoA and AoD computation for ue_ref with respect to current BS (no RT received, so it must be linked to old BS)
                    AoA = zeros(n_BS_all,1);
                    AoD = zeros(n_BS_all,1);
                    H = cell(n_BS_all,1);
                    H_params = cell(n_BS_all,1);
                    for bs_idx = 1:n_BS_all
                        pos_bs = BS_all(bs_idx,:);
                        [tmp_AoA tmp_AoD] = find_angles(ue_test_new_position,pos_bs);
                        AoA(bs_idx) = tmp_AoA; % AoA of 'ue' with respect to 'bs'
                        AoD(bs_idx) = tmp_AoD;
                        [H{bs_idx} H_params{bs_idx}] = compute_H_mobility(f,AoA(bs_idx),AoD(bs_idx),n_rx(antenna_idx),n_tx(antenna_idx)); % channel 
                    end
                end  
                
                if mod(t, T_tracking) < 1e-10 %% tracking is completed, select new beam pair --> moreover, all BSs possibly change their random orientations
                    
                    %% select new BS to connect to, after determining current PL (keep last computed PL map for the PL status, but with current distacne values)
                    BS_distance = abs( sqrt( (BS_all(:,1) - ue_test_new_position(1)).^2 + (BS_all(:,2) - ue_test_new_position(2)).^2 ) ); % new distance from BSs to current position of the UE
                    BS_PL = 10.^((32.4 + 21.6*log10(BS_distance) + 20*log10(f/1e9))/10)';
                    [serving_PL, serving_BS_idx] = min(BS_PL); % accoring to Assumption 3.4 and the definition of impariment PL in Eq. (2)
                    
                    %% AoA and AoD computation for ue_ref with respect to current BS (no RT received, so it must be linked to old BS)
                    AoA = zeros(n_BS_all,1);
                    AoD = zeros(n_BS_all,1);
                    H = cell(n_BS_all,1);                    
                    for bs_idx = 1:n_BS_all
                        pos_bs = BS_all(bs_idx,:);
                        [tmp_AoA tmp_AoD] = find_angles(ue_test_new_position,pos_bs);
                        AoA(bs_idx) = tmp_AoA; % AoA of 'ue' with respect to 'bs'
                        AoD(bs_idx) = tmp_AoD;
                        H{bs_idx} = compute_H_ssf(f,AoA(bs_idx),AoD(bs_idx),n_rx(antenna_idx),n_tx(antenna_idx), H_params{bs_idx}); % channel, keeping same small scale fading params
                    end
                    
                    %% All BSs (but the serving one) steer in random directions, while the UE is aligned with the serving BS
                    BF_vector_tx = cell(n_BS_all,1);
                    rnd_steering_angle = 2*pi*rand(n_BS_all,1); % random steering angles
                    for bs_idx = 1:n_BS_all
                        BF_vector_tx{bs_idx} = compute_BF_vector(n_tx(antenna_idx),rnd_steering_angle(bs_idx));
                    end
                    BF_vector_tx{serving_BS_idx} = compute_BF_vector(n_tx(antenna_idx),AoA(serving_BS_idx));
                    
                    BF_vector_rx = compute_BF_vector(n_rx(antenna_idx),AoD(serving_BS_idx));
                    Gain_num_new(index_internal) = (abs(conj(BF_vector_rx)*H{serving_BS_idx}*BF_vector_tx{serving_BS_idx}.').^2); % this is a test
                    
                    Gain_interference = zeros(n_BS_all,1);
                    for bs_idx = 1:n_BS_all
                        Gain_interference(bs_idx) = (abs(conj(BF_vector_rx)*H{bs_idx}*BF_vector_tx{bs_idx}.').^2);
                    end
                    
                    SINR_interference = Gain_interference( setdiff(1:n_BS_all , serving_BS_idx),: ) ./ ...
                        (BS_PL(setdiff(1:n_BS_all , serving_BS_idx))').*0.5;
                    
                    SINR_num =    Gain_num_new(index_internal)./ BS_PL(serving_BS_idx) ; % numerator of SINR (depends on the beamwidth)
                    SINR_den = sum(SINR_interference) + thermal_noise;
                    
                    avg_SINR(index_internal) = SINR_num ./ SINR_den;
                    rate(index_internal) = BW * log2(1+avg_SINR(index_internal));
                    index_internal = index_internal + 1;
                    
                else % stick with previous beam configuration, eventually change PL status. All BSs steer in the same directions as before
                    
                    BS_distance = abs( sqrt( (BS_all(:,1) - ue_test_new_position(1)).^2 + (BS_all(:,2) - ue_test_new_position(2)).^2 ) ); % new distance from BSs to current position of the UE
                    BS_PL = 10.^((32.4 + 21.6*log10(BS_distance) + 20*log10(f/1e9))/10)'; 
                    
                    %% AoA and AoD computation for ue_ref with respect to current BS (no RT received, so it must be linked to old BS)
                    AoA = zeros(n_BS_all,1);
                    AoD = zeros(n_BS_all,1);
                    H = cell(n_BS_all,1);                    
                    for bs_idx = 1:n_BS_all
                        pos_bs = BS_all(bs_idx,:);
                        [tmp_AoA tmp_AoD] = find_angles(ue_test_new_position,pos_bs);
                        AoA(bs_idx) = tmp_AoA; % AoA of 'ue' with respect to 'bs'
                        AoD(bs_idx) = tmp_AoD;
                        H{bs_idx} = compute_H_ssf(f,AoA(bs_idx),AoD(bs_idx),n_rx(antenna_idx),n_tx(antenna_idx), H_params{bs_idx}); % channel, keeping same small scale fading params
                    end
                       
                    %% DO NOT CHANGE BF VECTORS, but H matrix has changed
                    Gain_num_new(index_internal) = (abs(conj(BF_vector_rx)*H{serving_BS_idx}*BF_vector_tx{serving_BS_idx}.').^2); % this is a test

                    Gain_interference = zeros(n_BS_all,1);
                    for bs_idx = 1:n_BS_all
                        Gain_interference(bs_idx) = (abs(conj(BF_vector_rx)*H{bs_idx}*BF_vector_tx{bs_idx}.').^2);
                    end
                    
                    SINR_interference = Gain_interference( setdiff(1:n_BS_all , serving_BS_idx),: ) ./ ...
                        (BS_PL(setdiff(1:n_BS_all , serving_BS_idx))').*0.5;
                    
                    SINR_num =    Gain_num_new(index_internal)./ BS_PL(serving_BS_idx) ; % numerator of SINR (depends on the beamwidth)
                    SINR_den = sum(SINR_interference) + thermal_noise;
                    
                    avg_SINR(index_internal) = SINR_num ./ SINR_den;
                    rate(index_internal) = BW * log2(1+avg_SINR(index_internal));
                    index_internal = index_internal + 1;
                end
            end %% end simulation
            rate_final(n_iter) = mean(rate); % need to save for all values of vector_lambda, antennas
            
        end % end n iterations
    end % end vector_lambda
end % end antennas
stop_timer = toc;
fprintf('runtime: %4.3f s\n', stop_timer); 

formatSpec = '%8.5fGbps\n';
fileID = fopen('rate.txt','w');
fprintf(fileID,formatSpec,rate/1e9);