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

%% number of antennas for BS and UE
n_tx = 64; %WARNING: all values must be perfect squares
n_rx = 16; %WARNING: all values must be perfect squares

%% paramers for road
v = 100/3.6; % speed [m/s]
road_length = 1000;
d_R = 2.5; 
W_L = 3.7; % lane width
N_0 = 3; % # of obstacle lanes per direction
R = d_R + N_0*W_L; % road width
BS_per_km = 15;

%% parameters for simulation
n_rep_PL = 1;
theta_out = -5; %SINR outage threshold [dB]
outage_thresh = 10^(theta_out/10); %SINR outage threshold
T_sim = 1000/v; % simulation duration [s]
dt = 0.1; %  simulation step [s]
T_tracking = 0.1; % tracking periodicity for BF vector [s]
T_mul_users_update = 1; %how often BSs change n of connected users 
t_H = 0.3; % udpate of channel instances
n_users = 30; % mean number of users per BS, poisson r.v.

%% %%%%%%%%%%%%%%%%%%%%%%% Monte Carlo Method %%%%%%%%%%%%%%%%%%%%%%%%%%
rate_tmp = cell(n_rep_PL, 1);
tic; 
for iter = 1:n_rep_PL  

    fprintf('iteration: %d\n', iter);      
    
    ue_pos = [0, R - randi(3,1)*W_L + W_L/2, 0]; %(x,y,z) position of UE
    shared_data = BaseStation.sharedData;
    UE = UserEquipment(n_rx, f, ue_pos, v, T_tracking);
    shared_data.UE = UE;
    UE.sharedData = shared_data;
        
    n_BS_top = poissrnd(BS_per_km/1000 * (road_length),1,1); % number of BSs
    n_BS_bottom = poissrnd(BS_per_km/1000 * (road_length),1,1); % number of BSs   
    BS_distance_avg_top = 1000/n_BS_top;
    BS_distance_avg_bottom = 1000/n_BS_bottom;
    delta_top = ceil(BS_distance_avg_top/8);
    delta_bottom = ceil(BS_distance_avg_bottom/8);    
    
    allBS = cell(n_BS_top + n_BS_bottom, 1);
    for i = 1:n_BS_top
        pos = [(i-1)*BS_distance_avg_top+unifrnd(-delta_top , delta_top) , 2 * R, 8];
        allBS{i} = BaseStation(i, n_tx, f, BW, pos, t_H, T_tracking, n_users); 
    end
    
    for i = 1:n_BS_bottom
        pos = [(i-1)*BS_distance_avg_bottom+unifrnd(-delta_bottom , delta_bottom) , 0, 8];
        allBS{i + n_BS_top} = BaseStation(i + n_BS_top, n_tx, f, BW, pos, t_H, T_tracking, n_users); 
    end
    
%     %% for debug, plot BS disposition
%     figure;
%     hold on;
%     for i = 1:n_BS_top+n_BS_bottom
%         plot(allBS{i}.pos(1)+1, allBS{i}.pos(2), '*')
%         text(allBS{i}.pos(1)+1, allBS{i}.pos(2), int2str(allBS{i}.ID));
%     end
%     hold off;
%     %%
    
    shared_data.servingBS = allBS{1}; %just for initialization
    UE.init();
    for i = 1:length(allBS)
            allBS{i}.init();
    end
    index_internal = 1;
    rate = zeros(1, length(dt:dt:T_sim));
    for t = dt:dt:T_sim %start simulation
        UE.update(t, dt);
        
        for i = 1:length(allBS)
            allBS{i}.update(t, dt);
        end
        
%         %% for debug
%          fprintf('-\n');
%         fprintf('serving BS: %d\n', shared_data.servingBS.ID);
%         fprintf('UE xpos: %3.3f\n', UE.pos(1));
%         id = 1;
%         for i = 1:length(allBS)
%             if allBS{i}.signal_power_at_ue > allBS{id}.signal_power_at_ue
%                 id = i;
%             end
%         end
%         fprintf('max power BS id: %d\n', id);
%         %%
     
        SINR_num = shared_data.servingBS.signal_power_at_ue; % numerator of SINR (depends on the beamwidth)
        
        SINR_interference = -SINR_num; % we don't want such signal_power_at_ue here, but in the for loop we sum it for "error"
        for i = 1:length(allBS)
            SINR_interference = SINR_interference + (allBS{i}.signal_power_at_ue);    
        end       
        SINR_interference = SINR_interference * 0.5; % ??? why is it divided by 2 ???
                
        SINR_den = SINR_interference + thermal_noise;

        SINR = SINR_num / SINR_den; 
        if SINR < outage_thresh
            SINR = 0; %outage -> no connection between BS and UE
        end
        
%         %% for debug          
%         D(index_internal) = norm(UE.pos - shared_data.servingBS.pos) / 1000;
%         PL(index_internal) = shared_data.servingBS.PL;
%         G(index_internal) = shared_data.servingBS.signal_power_at_ue / shared_data.servingBS.PL;
%         ASINR(index_internal) = SINR;
%         %%
        
        %consider rate for the n-users loaded BS
        rate(index_internal) = shared_data.servingBS.BW * log2(1+SINR) / shared_data.servingBS.n;   %OUTPUT of this Monte Carlo iteration
        index_internal = index_internal + 1; 
    end
    
%     %% for debug
%     figure;
%     plot(PL);
%     title('Path Loss');
%     figure;
%     plot(G);
%     title('signal power at ue');
%     figure;
%     plot(D);
%     title('Distance');   
%     figure;
%     plot(ASINR);
%     title('AVG SINR');
%     figure;
%     plot(rate);
%     title('Rate');
%     %%
        
    rate_tmp{iter} = [min(rate), mean(rate), max(rate), std(rate)]; % entire OUTPUT of the Monte Carlo method
end

stop_timer = toc; 
fprintf('runtime: %4.3f s\n', stop_timer);  

%% %%%%%%%%%%%%%%%%%%%%%%reorder output%%%%%%%%%%%%%%%%%%%%%%
min_rate_final = zeros(1,n_rep_PL);
mean_rate_final = zeros(1,n_rep_PL); 
max_rate_final = zeros(1,n_rep_PL);
std_rate_final = zeros(1,n_rep_PL);

for iter =  1:n_rep_PL 
    min_rate_final(1, iter) = rate_tmp{iter}(1);
    mean_rate_final(1, iter) = rate_tmp{iter}(2);
    max_rate_final(1, iter) = rate_tmp{iter}(3);
    std_rate_final(1, iter) = rate_tmp{iter}(4);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

formatSpec = '%4d->\tmin\t%2.6f\tmean\t%2.6f\tmax\t%2.6f\tstd\t%2.6f\tGbps\n';  %tab separeted values
tmp = 1:n_rep_PL;
iter_vec = 1:n_rep_PL;

tofile =    [   iter_vec;    ...
                min_rate_final / 1e9;  ...
                mean_rate_final / 1e9;  ...
                max_rate_final / 1e9;  ...
                std_rate_final / 1e9   ...
            ];

fname = sprintf('rate_%d_%d_%d.txt', BS_per_km, n_tx, n_rx);
fileID = fopen(fname,'w');
fprintf(fileID, 'configuration:\t%d BSperKm\t%d n_tx\t%d n_rx\n', [BS_per_km, n_tx, n_rx]);
fprintf(fileID, formatSpec, tofile);
fclose(fileID);

%print a fast report of the simulation
fprintf('mean rate: %2.6f Gbps\n', mean(mean_rate_final(1, :)) / 1e9);