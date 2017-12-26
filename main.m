%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%Simulate a mmWave scenario with proactive caching.
%
%Authors: Frizziero - Suman - Dell'Eva
%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all;
close all;
clc;

seed = sum(clock);
rng(seed);

print_scenario = false;
addpath(genpath('fnc'))
addpath(genpath('utils'))
addpath(genpath('Classes'))

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
    n_tx = n_tx_array(antenna_idx);
    n_rx = n_rx_array(antenna_idx);
    fprintf('iteration: %d [%d/%d -> %d/%d]\n', iter, antenna_idx, length(n_tx_array),n_iter, n_rep_PL);
    
    ue_pos = [rand(1,1) * (road_length - usefull_road_length), R - randi(3,1)*W_L + W_L/2];
    shared_data = BaseStation.sharedData;
    UE = UserEquipment(n_rx, ue_pos, v, T_tracking);
    shared_data.UE = UE;
    UE.sharedData = shared_data;
    
    
    n_BS_top = poissrnd(vector_lambda_bs * (road_length),1,1); % number of BSs
    BS_top = [unifrnd(0, road_length, n_BS_top, 1) , 2 * R * ones(n_BS_top,1)]; %position of BSs

    n_BS_bottom = poissrnd(vector_lambda_bs * (road_length),1,1); % number of BSs
    BS_bottom = [unifrnd(0, road_length, n_BS_bottom, 1) , zeros(n_BS_bottom,1)]; %position of BSs
             
    allBS = cell(n_BS_top + n_BS_bottom, 1);
    for i = 1:n_BS_top
        allBS{i} = BaseStation(n_tx, f, BW, BS_top(i, :), t_H, T_tracking, n_users); 
    end
    
    for i = 1:n_BS_bottom
        allBS{i + n_BS_top} = BaseStation(n_tx, f, BW, BS_bottom(i, :), t_H, T_tracking, n_users); 
    end
    
    shared_data.servingBS = allBS{1}; %just for initialization
    UE.init();
    for i = 1:length(allBS)
            allBS{i}.init();
    end
    index_internal = 1;
    rate = zeros(1, length(0.1:t_offset:T_sim));
    for t = t_offset:t_offset:T_sim %start simulation
        UE.update(t);
        
        for i = 1:length(allBS)
            allBS{i}.update(t);
        end
    
    
        SINR_interference = -(shared_data.servingBS.GAIN / shared_data.servingBS.PL) * 0.5; % we wont such gain here, but in the for loop we sum it for "error"
        for i = 1:length(allBS)
            SINR_interference = SINR_interference + (allBS{i}.GAIN / allBS{i}.PL) * 0.5;    
        end

        SINR_num = shared_data.servingBS.GAIN / shared_data.servingBS.PL ; % numerator of SINR (depends on the beamwidth)
        SINR_den = SINR_interference + thermal_noise;

        avg_SINR = SINR_num ./ SINR_den; 
        
        %consider rate for the n-users loaded BS
        rate(index_internal) = shared_data.servingBS.BW * log2(1+avg_SINR) / shared_data.servingBS.n;   %OUTPUT of this Monte Carlo iteration
        index_internal = index_internal + 1; 
    end
    
    rate_tmp{iter} = [min(rate), mean(rate), max(rate), std(rate)]; % entire OUTPUT of the Monte Carlo method
end

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
                [min_rate_final(1, :), min_rate_final(2, :), min_rate_final(3, :) ]  / 1e9;  ...
                [mean_rate_final(1, :), mean_rate_final(2, :), mean_rate_final(3, :) ]  / 1e9;  ...
                [max_rate_final(1, :), max_rate_final(2, :), max_rate_final(3, :) ]  / 1e9;  ...
                [std_rate_final(1, :), std_rate_final(2, :), std_rate_final(3, :) ]  / 1e9   ...
            ];

fname = sprintf('rate_%dBSperKM.txt', BS_per_km);
fileID = fopen(fname,'w');
fprintf(fileID, formatSpec, tofile);
fclose(fileID);

%print a fast report of the simulation
fprintf('antenna configuration 1 mean rate: %2.6f Gbps\n', mean(mean_rate_final(1, :)) / 1e9);
fprintf('antenna configuration 2 mean rate: %2.6f Gbps\n', mean(mean_rate_final(2, :)) / 1e9);
fprintf('antenna configuration 3 mean rate: %2.6f Gbps\n', mean(mean_rate_final(3, :)) / 1e9);