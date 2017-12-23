%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% V2I in vehicular network operating at millimeter wave.
% A multiu-lane highway scenario is considered.
%
% Marco Giordani
% July 2017
%revised by Frizziero - Suman - Dell'Eva
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all;
%close all;
clc;
rng(4)

print_scenario = false;
addpath(genpath('fnc'))
addpath(genpath('utils'))

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

n_tx_array = [4 16 64]; %WARNING: all values must be perfect squares
n_rx_array = [4  4 16]; %WARNING: all values must be perfect squares


%% parameters for simulation
n_rep_PL = 10;
theta_out = -5; %SINR outage threshold [dB]
theta_out_lin = 10.^(theta_out./10); %SINR outage threshold
T_sim = 20; % simulation duration [s]
t_offset = 0.1; %  simulation step [s]
T_tracking = 0.1; % tracking periodicity [s]
T_mul_users_update = 1; %how often BSs change n of connected users 
t_H = 0.3; % udpate of channel instances
n_users = 30; % mean number of users per BS, poisson r.v.

%% paramers for road
v = 100/3.6; % speed [m/s]
road_length = 1000;
usefull_road_length = v * T_sim; % positive length of the road [m] %modified by Frizziero
road_start = -1000; % negative length of the road [m]
d_R = 2.5; 
W_L = 3.5; % truck length
N_0 = 3; % # of obstacle lanes per direction
R = d_R + N_0*W_L; % road width
BS_per_km = 15;

%% %%%%%%%%%%%%%%%%%%%%%%% set up parallelization %%%%%%%%%%%%%%%%%%%%%%%%%%
permutation = randperm(n_rep_PL*length(n_tx_array)); %to equally spread workload on workers
vector_lambda_bs = BS_per_km/1000;
tmp = 1:n_rep_PL;
iter_vec = [];
for i = 1:length(n_tx_array)
    iter_vec = [iter_vec, tmp];
end
antenna_idx_vec = [];
tmp = ones(1, n_rep_PL);
for i = 1:length(n_tx_array)
    antenna_idx_vec = [antenna_idx_vec, tmp];
    tmp = tmp + 1;
end

iteration_map = (1:n_rep_PL)';
for i = 1:length(n_tx_array)-1
    iteration_map = [iteration_map; (1:n_rep_PL)'];
    tmp = tmp + 1;
end

antenna_map = ones(n_rep_PL, 1);
for i = 1:length(n_tx_array)-1
    antenna_map = [antenna_map; ones(n_rep_PL, 1) + i];
end
iteration_map = iteration_map(permutation);
antenna_map = antenna_map(permutation);
rate_tmp = cell(length(n_tx_array) * n_rep_PL, 1);
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

try
    nc = feature('numcores');
    if nc > 32
        nc = 32; %max allowed cpu cores usage for DEI's cluster
    end
	parpool(nc)
catch E
    %parpool is already up and running
end

tic; 
parfor iter = 1:n_rep_PL*length(n_tx_array)  

    n_iter = iteration_map(iter);    
    antenna_idx = antenna_map(iter);
    n_tx = n_tx_array(antenna_idx)
    n_rx = n_rx_array(antenna_idx)
    fprintf('iteration: %d [%d/%d -> %d/%d]\n', iter, antenna_idx, length(n_tx_array),n_iter, n_rep_PL);

    UE = [rand(1,1) * (road_length - usefull_road_length), R - randi(3,1)*W_L + W_L/2]; % choose one of the three lanes randomly %added by Frizziero

    %% deployment TOP part of road
    obs_lanes_top = cell(1,N_0);

    n_BS_top = poissrnd(vector_lambda_bs * (road_length),1,1); % number of BSs
    BS_top = [unifrnd(0, road_length, n_BS_top,1) , 2*R * ones(n_BS_top,1)]; %position of BSs

    n_BS_bottom = poissrnd(vector_lambda_bs * (road_length),1,1); % number of BSs
    BS_bottom = [unifrnd(0, road_length, n_BS_bottom,1) , 0 * ones(n_BS_bottom,1)]; %position of BSs

    BS_all = [BS_top; BS_bottom];

    n_BS_all = n_BS_top + n_BS_bottom
    BS_distance = abs( sqrt( (BS_all(:,1) - UE(1)).^2 + (BS_all(:,2) - UE(2)).^2 ) ); % new distance from BSs to current position of the UE

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% PRINT SCENARIO (if enabled)
%     if print_scenario
%         %hold on
%         scatter(BS_all(:,1),BS_all(:,2),'*r'); hold on;
%         scatter(UE(1),UE(2),'dm')
%     end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %% Find PL and PL states, according to the distance
    BS_PL = 10.^((32.4 + 21.6*log10(BS_distance) + 20*log10(f/1e9))/10)';

    %% MAKE ATTACHMENT DECISION (max. imparing PL, accoridng to the PL defined in Eq. (2))
    [serving_PL, serving_BS_idx] = min(BS_PL); % accoring to Assumption 3.4 and the definition of impariment PL in Eq. (2)
    distance_serving_BS = abs(BS_distance(serving_BS_idx)); % distance from UE in (0,0) to current serving BS

    %% AoA, AoD, spatial signatures, H computation for ues with respect to all BSs in the system
    H = cell(n_BS_all,1);
    H_params = cell(n_BS_all,1);
    for bs_idx = 1:n_BS_all
        pos_bs = BS_all(bs_idx,:);
        [AoA, AoD] = find_angles(UE,pos_bs);
        [H{bs_idx}, H_params{bs_idx}] = compute_H_mobility(f, AoA, AoD, n_rx, n_tx, BW); % channel 
    end

    %% All BSs (but the serving one) steer in random directions, while the UE is aligned with the serving BS
    pos_bs = BS_all(serving_BS_idx,:);
    [AoA, AoD] = find_angles(UE,pos_bs);
    [BF_vector_tx, BF_vector_rx] = update_BF_vectors(AoA, AoD, n_tx, n_rx, n_BS_all, serving_BS_idx);
                                                            
    Gain_num = (abs(conj(BF_vector_rx)*H{serving_BS_idx}*BF_vector_tx{serving_BS_idx}.').^2); 

    Gain_interference = zeros(n_BS_all,1);
    for bs_idx = 1:n_BS_all
        Gain_interference(bs_idx) = (abs(conj(BF_vector_rx)*H{bs_idx}*BF_vector_tx{bs_idx}.').^2);
    end

    SINR_interference = Gain_interference( setdiff(1:n_BS_all , serving_BS_idx),: ) ./ ...
        (BS_PL(setdiff(1:n_BS_all , serving_BS_idx))').*0.5;

    SINR_num =    Gain_num ./ BS_PL(serving_BS_idx)  ; % numerator of SINR (depends on the beamwidth)

    SINR_den = sum(SINR_interference) + thermal_noise;
    BW_shared = BW / poissrnd(n_users); % simulate load on BS due to multiple concurrent users
    %original_avg_SINR = SINR_num ./ SINR_den; % ??? never used after or before
    %original_rate = BW * log2(1+original_avg_SINR); % ??? never used after or before

    %% consider time-varying simulation
    ue_test_new_position = [UE(1) UE(2)];
    rate = zeros(1, length(0.1:t_offset:T_sim)) 
    index_internal = 1; %internal cycle index
    for t = t_offset:t_offset:T_sim %t_offset:t_offset:T_sim

        space_offset = t_offset * v;
        ue_test_new_position = [ue_test_new_position(1)+space_offset ue_test_new_position(2)]; %new position of the user, after it has moved

        if mod(t, t_H) < 1e-10 %% change PL status for all BSs     
            %% AoA and AoD computation for ue_ref with respect to current BS (no RT received, so it must be linked to old BS)            
            H = cell(n_BS_all,1);
            H_params = cell(n_BS_all,1);
            for bs_idx = 1:n_BS_all
                pos_bs = BS_all(bs_idx,:);
                [AoA, AoD] = find_angles(ue_test_new_position, pos_bs);
                [H{bs_idx}, H_params{bs_idx}] = compute_H_mobility(f, AoA, AoD, n_rx, n_tx, BW); % channel 
            end
        end  

        if mod(t, T_tracking) < 1e-10 %% tracking is completed, select new beam pair --> moreover, all BSs possibly change their random orientations
            %% select new BS to connect to, after determining current PL (keep last computed PL map for the PL status, but with current distacne values)
            BS_distance = abs( sqrt( (BS_all(:,1) - ue_test_new_position(1)).^2 + (BS_all(:,2) - ue_test_new_position(2)).^2 ) ); % new distance from BSs to current position of the UE
            BS_PL = 10.^((32.4 + 21.6*log10(BS_distance) + 20*log10(f/1e9))/10)';
            [serving_PL, serving_BS_idx] = min(BS_PL); % accoring to Assumption 3.4 and the definition of impariment PL in Eq. (2)

            %% AoA and AoD computation for ue_ref with respect to current BS (no RT received, so it must be linked to old BS)           
            H = cell(n_BS_all,1);                    
            for bs_idx = 1:n_BS_all
                pos_bs = BS_all(bs_idx,:);
                [AoA, AoD] = find_angles(ue_test_new_position, pos_bs);
                H{bs_idx} = compute_H_ssf(f, AoA, AoD, n_rx, n_tx, H_params{bs_idx}, BW); % channel, keeping same small scale fading params
            end

            %% All BSs (but the serving one) steer in random directions, while the UE is aligned with the serving BS
            pos_bs = BS_all(serving_BS_idx,:);
            [AoA, AoD] = find_angles(ue_test_new_position, pos_bs);
            [BF_vector_tx, BF_vector_rx] = update_BF_vectors(AoA, AoD, n_tx, n_rx, n_BS_all, serving_BS_idx);                

        else % stick with previous beam configuration, eventually change PL status. All BSs steer in the same directions as before
            %% DO NOT CHANGE BF VECTORS, but H matrix has changed 

            %% AoA and AoD computation for ue_ref with respect to current BS (no RT received, so it must be linked to old BS)
            H = cell(n_BS_all,1);                    
            for bs_idx = 1:n_BS_all
                pos_bs = BS_all(bs_idx,:);
                [AoA, AoD] = find_angles(ue_test_new_position,pos_bs);
                H{bs_idx} = compute_H_ssf(f, AoA, AoD, n_rx, n_tx, H_params{bs_idx}, BW); % channel, keeping same small scale fading params
            end                   
        end
        
        Gain_num_new = (abs(conj(BF_vector_rx)*H{serving_BS_idx}*BF_vector_tx{serving_BS_idx}.').^2); % this is a test %modified by Frizziero

        Gain_interference = zeros(n_BS_all,1);
        for bs_idx = 1:n_BS_all
            Gain_interference(bs_idx) = (abs(conj(BF_vector_rx)*H{bs_idx}*BF_vector_tx{bs_idx}.').^2);
        end

        SINR_interference = Gain_interference( setdiff(1:n_BS_all , serving_BS_idx),: ) ./ ...
            (BS_PL(setdiff(1:n_BS_all , serving_BS_idx))').*0.5;

        SINR_num = Gain_num_new ./ BS_PL(serving_BS_idx) ; % numerator of SINR (depends on the beamwidth)
        SINR_den = sum(SINR_interference) + thermal_noise;
        
        avg_SINR = SINR_num ./ SINR_den; 
        if mod(t, T_mul_users_update) < 1e-10 %update every 1 sec
            BW_shared = BW / poissrnd(n_users); % simulate load on BS due to multiple concurrent users
        end
        rate(index_internal) = BW_shared * log2(1+avg_SINR);   %OUTPUT of this Monte Carlo iteration
        index_internal = index_internal + 1;    
    end %% end simulation

    rate_tmp{iter} = [min(rate), mean(rate), max(rate), std(rate)]; % entire OUTPUT of the Monte Carlo method
end % end n iterations

stop_timer = toc; 
fprintf('runtime: %4.3f s\n', stop_timer);  

%% %%%%%%%%%%%%%%%%%%%%%%reorder output%%%%%%%%%%%%%%%%%%%%%%
min_rate_final = zeros(length(n_tx_array),n_rep_PL);
mean_rate_final = zeros(length(n_tx_array),n_rep_PL); 
max_rate_final = zeros(length(n_tx_array),n_rep_PL);
std_rate_final = zeros(length(n_tx_array),n_rep_PL);

for iter =  1:n_rep_PL*length(n_tx_array)
    antenna_idx = antenna_map(iter);
    n_iter = iteration_map(iter);
    min_rate_final(antenna_idx, n_iter) = rate_tmp{iter}(1);
    mean_rate_final(antenna_idx, n_iter) = rate_tmp{iter}(2);
    max_rate_final(antenna_idx, n_iter) = rate_tmp{iter}(3);
    std_rate_final(antenna_idx, n_iter) = rate_tmp{iter}(4);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

formatSpec = '%d|%d->\tmin\t%2.6f\tmean\t%2.6f\tmax\t%2.6f\tstd\t%2.6f\tGbps\n';  %tab separeted values
tmp = 1:n_rep_PL;
iter_vec = [];
for i = 1:length(n_tx_array)
    iter_vec = [iter_vec, tmp];
end
antenna_idx_vec = [];
tmp = ones(1, n_rep_PL);
for i = 1:length(n_tx_array)
    antenna_idx_vec = [antenna_idx_vec, tmp];
    tmp = tmp + 1;
end

tofile =    [   antenna_idx_vec; iter_vec;    ...
                reshape(min_rate_final / 1e9, [1, length(iter_vec)]);  ...
                reshape(mean_rate_final / 1e9, [1, length(iter_vec)]);  ...
                reshape(max_rate_final / 1e9, [1, length(iter_vec)]);  ...
                reshape(std_rate_final / 1e9, [1, length(iter_vec)])   ...
            ];

fname = sprintf('rate_%dBSperKM.txt', BS_per_km);
fileID = fopen(fname,'w');
fprintf(fileID, formatSpec, tofile);