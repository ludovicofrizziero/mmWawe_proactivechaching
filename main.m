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
%v = 100/3.6; % speed [m/s]
road_length = 1000;
d_R = 2.5; 
W_L = 3.75; % lane width
N_0 = 3; % # of obstacle lanes per direction
R = d_R + N_0*W_L; % road width
BS_per_km = 5; %BSs' *expected* density. Just for one row (top/bottom) -> actual density is 2*BS_per_km

%% parameters for simulation
n_rep_PL = 300; % number of repetition per choice of parameters
theta_out = -5; %SINR outage threshold [dB]
outage_thresh = 10^(theta_out/10); %SINR outage threshold
% T_sim = floor(1000/v); % simulation duration [s]
dt = 0.1; %  simulation step [s]
T_tracking = 0.1; % tracking periodicity for BF vector [s]
T_mul_users_update = 1; %how often BSs change n of connected users 
t_H = 0.3; % udpate of channel instances
n_users = 25; % mean number of users per BS, poisson r.v.
alloc_func = {@custom_solver3; @random_allocation3; @random_allocation1}


DEBUG = n_rep_PL < 2;
SAVE_DATA_VERBOSE = true;

%% deploy BS
[n_BS_top, n_BS_bottom, pos_top, pos_bottom] = deploy_bs(BS_per_km, road_length, R);
%%

fprintf('Going to delete all previous results and start a new batch of simulations. PRESS ANY KEY TO CONTINUE.\n')
if ~DEBUG
    pause();
end
delete('RESULTS//savings_v*');
global_start = tic;
for v = (70:10:130)/3.6 %set of velocities for the ue [m/s]  
%     [n_BS_top, n_BS_bottom, pos_top, pos_bottom] = deploy_bs(BS_per_km, road_length, R);
%     for ue_r = 0.068 %[0.016, 0.045, 0.068] %set of rates for the ue [Gbit/s]        
        for alloc_func_idx = 1:max(size(alloc_func))
            fprintf('v = %3.0f Km/h, alloc_func = %s\n', round(v*3.6), func2str(alloc_func{alloc_func_idx}));
            T_sim = floor(1000/v); % simulation duration [s]

            %% %%%%%%%%%%%%%%%%%%%%%%% Monte Carlo Method %%%%%%%%%%%%%%%%%%%%%%%%%%             
            rate_tmp = cell(n_rep_PL, 1);        
            savings = cell(n_rep_PL, 1);
            tic; 
            parfor iter = 1:n_rep_PL                                       
                ue_r = (0.15 - 0.06) * rand(1) + 0.06; %Gbit/s
                start_ue_pos = [0, d_R + randi(3,1)*W_L - W_L/2, 0]; %(x,y,z) position of UE
                shared_data = BaseStation.sharedData;
                UE = UserEquipment(n_rx, f, start_ue_pos, v, ue_r, T_tracking, T_sim/dt);
                shared_data.UE = UE;
                UE.sharedData = shared_data; 
                fprintf('\titeration: %d\n' , iter);

                ok = false;
                allBS = [];
                X = [];
                chunks = [];
                n_retry = 0;
                while ~ok && n_retry < 100
                    [n_BS_top, n_BS_bottom, pos_top, pos_bottom] = deploy_bs(BS_per_km, road_length, R);

                    allBS = cell(n_BS_top + n_BS_bottom, 1);
                    for i = 1:n_BS_top                
                        allBS{i} = BaseStation(2*i, n_tx, f, BW, pos_top(i, :), t_H, T_tracking, n_users); 
                    end

                    for i = 1:n_BS_bottom
                        allBS{i + n_BS_top} = BaseStation(2*i-1, n_tx, f, BW, pos_bottom(i, :), t_H, T_tracking, n_users); 
                    end    

                    shared_data.servingBS = allBS{1}; %just for initialization
                    UE.init();
                    for i = 1:length(allBS)
                            allBS{i}.init();
                    end

                    %% allocate file to BS
                    [X, chunks, ok] = alloc_func{alloc_func_idx}(allBS, UE, BS_per_km, DEBUG); %max(n_BS_bottom, n_BS_top), DEBUG);
                    %%     
                    n_retry = n_retry + 1;
                end

                for i = 1:size(X)
                    allBS{i}.allocate_memory_for_ue(chunks(i) * X(i));
                end

                if DEBUG
                    disp(int8(X'));
                end
                                        

                %% start simulation
                index_internal = 1;
                sim_steps = dt:dt:T_sim;
                rate = zeros(length(sim_steps), 1);
                servingBS_IDs = zeros(length(sim_steps), 1);
                for t = sim_steps 
                    distance_handover(allBS, shared_data); 
%                     power_handover(allBS, shared_data); 
                    UE.update(t, dt);

                    for i = 1:length(allBS)
                        allBS{i}.update(t, dt);
                    end

                    SINR_num = shared_data.servingBS.signal_power_at_ue; % numerator of SINR (depends on the beamwidth)

                    SINR_interference = -SINR_num; % we don't want such signal_power_at_ue here, but in the for loop we sum it for "error"
                    for i = 1:length(allBS)
                        SINR_interference = SINR_interference + (allBS{i}.signal_power_at_ue);    
                    end       
                    SINR_interference = SINR_interference * 0.5; % ??? why is it divided by 2 ???

                    SINR_den = SINR_interference + thermal_noise/shared_data.servingBS.n; %we must account for the avaiable bandwith to the UE, not the total BW

                    SINR = SINR_num / SINR_den; 
                    if SINR < outage_thresh
                        SINR = 0; %outage -> no connection between BS and UE
                    end       

                    %% file transmission
                    %consider rate for the n-users loaded BS
                    overhead_performance_loss = 1 - 0.3; %consider 30% of performance loss due to modulation / packet headers ...
                    r = shared_data.servingBS.BW * log2(1+SINR) / shared_data.servingBS.n;% * overhead_performance_loss;       

                    f_chunk = shared_data.servingBS.download_file(dt, r);
                    UE.receive_file_chunk(f_chunk);
                    %%

                    rate(index_internal) = r;
                    servingBS_IDs(index_internal) = shared_data.servingBS.ID;
                    index_internal = index_internal + 1; 
                end    

                if SAVE_DATA_VERBOSE
                    N = max(size(allBS));
                    [ue_buffer, ue_max_buffer, ue_lost_data, ue_waiting_time] = UE.dump_data();
                    BSs_pos = zeros(N, 3);
                    BSs_mem_state = zeros(N, 1);
                    all_ids = zeros(N,1);
                    for i = 1:N
                        BSs_pos(i, :) = allBS{i}.pos;
                        BSs_mem_state(i) = allBS{i}.memory;
                        all_ids(i) = allBS{i}.ID;
                    end        

                    s = struct;
                    s.ue_buffer = ue_buffer;
                    s.ue_max_buffer = ue_max_buffer;
                    s.ue_lost_data = ue_lost_data;
                    s.ue_waiting_time = ue_waiting_time;
                    s.ue_requested_rate = UE.requested_rate;
                    s.ue_velocity = UE.m_vel;
                    s.rate = rate;
                    s.servingBS_IDs = servingBS_IDs;
                    s.BSs_pos = BSs_pos;
                    s.BSs_mem_state = BSs_mem_state;
                    s.all_ids = all_ids;
                    s.X = X;
                    s.chunks = chunks .* X;        

                    savings{iter} = s;
                end

                rate_tmp{iter} = [min(rate), mean(rate), max(rate), std(rate)]; % entire OUTPUT of the Monte Carlo method
            end

            stop_timer = toc; 
            fprintf('\truntime: %4.3f s\n', stop_timer);  

%             %% %%%%%%%%%%%%%%%%%%%%%%reorder output%%%%%%%%%%%%%%%%%%%%%%
%             min_rate_final = zeros(1,n_rep_PL);
            mean_rate_final = zeros(1,n_rep_PL); 
%             max_rate_final = zeros(1,n_rep_PL);
%             std_rate_final = zeros(1,n_rep_PL);
% 
            for iter =  1:n_rep_PL 
%                 min_rate_final(1, iter) = rate_tmp{iter}(1);
                mean_rate_final(1, iter) = rate_tmp{iter}(2);
%                 max_rate_final(1, iter) = rate_tmp{iter}(3);
%                 std_rate_final(1, iter) = rate_tmp{iter}(4);
            end
%             %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
%             formatSpec = '%4d->\tmin\t%2.6f\tmean\t%2.6f\tmax\t%2.6f\tstd\t%2.6f\tGbps\n';  %tab separated values
%             tmp = 1:n_rep_PL;
%             iter_vec = 1:n_rep_PL;
% 
%             tofile =    [   iter_vec;    ...
%                             min_rate_final / 1e9;  ...
%                             mean_rate_final / 1e9;  ...
%                             max_rate_final / 1e9;  ...
%                             std_rate_final / 1e9   ...
%                         ];
% 
%             fname = sprintf('RESULTS\\rate_%d_%d_%d.txt', BS_per_km, n_tx, n_rx);
%             fileID = fopen(fname,'w');
%             fprintf(fileID, 'configuration:\t%d BSperKm\t%d n_tx\t%d n_rx\n', [BS_per_km, n_tx, n_rx]);
%             fprintf(fileID, formatSpec, tofile);
%             fclose(fileID);

            if SAVE_DATA_VERBOSE
                vel = num2str(round(v*3.6));
                if v*3.6 < 100
                    vel = strcat('0', vel);
                end
                name = strcat('RESULTS//savings_v', vel, '_func', num2str(alloc_func_idx), '.mat');
                save(   name, 'savings' );
            end
            
            %print a fast report of the simulation
            fprintf('\tmean rate: %2.6f Gbps\n', mean(mean_rate_final(1, :)) / 1e9);
        end
%     end
end
global_stop = toc(global_start);
fprintf('total runtime: %6.3f min\n', global_stop/60);
