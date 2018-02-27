% clear all;
% close all;
% clc;
% 
% all_files = dir('RESULTS//savings*');
% load(strcat('RESULTS//', all_files(1).name));
% 
% s = savings{1};
% 
% figure;
% plot(s.servingBS_IDs, '-');
% title('serving BS ID');
% 
% figure;
% hold on
% plot(s.all_ids, s.chunks, 'o');
% plot(s.all_ids, s.BSs_mem_state, '*');
% title('memory state at BSs');
% legend('starting situation', 'final situation');
% hold off
% 
% figure;
% plot(s.rate, '-');
% title('rate [bit/s]');
% 
% figure;
% hold on
% plot(s.ue_buffer);
% plot([1,max(size(s.ue_buffer))], ones(2,1) * double(s.ue_max_buffer));
% hold off
% title(strcat('UE buffer (with consumption rate of', num2str(s.ue_requested_rate, ' %1.3f'), ' Gbps )'));
% 
% figure;
% plot(s.ue_waiting_time, 'o');
% title('UE waiting times');
% 
% figure;
% plot(s.ue_lost_data);
% title('UE lost data');
% fprintf('PAUSED, PRESS ANY KEY TO CONTINUE\n');
% pause();

close all;
clear all;
addpath(genpath('utils'))

out = cell(1,1);
vels = 70:10:130;
for v = vels
    vel = num2str(v);
    if v < 100
        vel = strcat('0', vel);
    end
    name = strcat('RESULTS//savings_v', vel, '*');
    all_files = dir(name);
    for func_idx = 1:max(size(all_files))      
        saves = load(strcat('RESULTS//', all_files(func_idx).name));
        saves = saves.savings; %due to load() strangeness
        ue_d_lost = [];
        ue_w_time = [];
        bs_mem_left = [];        
        tot_chunks = [];
        ue_rates = [];
        for i = 1:max(size(saves)) % # of monte carlo iterations per file
            s = saves{i};
            ue_d_lost = [ue_d_lost; sum(s.ue_lost_data(1:end-5))];%the last UE's positions are degenerate due to simulation ending  
            ue_w_time = [ue_w_time; sum(s.ue_waiting_time(1:end-5))];
            bs_mem_left = [bs_mem_left; sum(s.BSs_mem_state(1:end-15))]; 
            tot_chunks = [tot_chunks; sum(s.chunks)];
            ue_rates = [ue_rates ; s.ue_requested_rate];
            if sum(s.chunks) == 0
                tot_chunks(end) = 1e9; %sometimes a solution is not found since all BS have full memory
            end
        end
        
        idx = (v-70)/10 + 1;
        out{idx, func_idx} = struct;
        out{idx, func_idx}.mean_lost = mean(ue_d_lost);        
        out{idx, func_idx}.mean_time = mean(ue_w_time);
        out{idx, func_idx}.CI_time = ConfIntervals(ue_w_time);
        out{idx, func_idx}.mean_mem = mean(bs_mem_left);
        out{idx, func_idx}.time_verbose = ue_w_time;
        out{idx, func_idx}.lost_verbose = ue_d_lost + bs_mem_left;
        out{idx, func_idx}.chunks_verbose = tot_chunks;
        out{idx, func_idx}.sim_duration = length(s.ue_waiting_time) * 0.1;
        out{idx, func_idx}.ue_rates = ue_rates;
        %out{idx, func_idx}.ue_buff = saves{1}.ue_buffer;
        %out{idx, func_idx}.ue_max_buff = saves{1}.ue_max_buffer;
    end
end

symbols = {'b-*', 'r-^', 'k-o', 'g-x'}; %for functions
colors = {'b:', 'r:', 'k:', 'g:'}; %for Conf. Interv.

%% lost memory
figure;
title('Average lost data due to buffer overflow');
hold on
grid on
for i = 1:min(size(out))
    y = [];
    for j = 1:max(size(out))
        y = [y; out{j, i}.mean_lost];
    end
    plot(vels, y / 8e6, symbols{i});
end
hold off
xlabel('velocity [km/h]');
ylabel('[MBytes]');
legend('Custom', 'Random 1', 'Random 2');
%%

%% wait time
figure;
title('Average cumulative wait time')
hold on
grid on
legend_subset = [];
for i = 1:min(size(out))
    y = [];
    ci = [];
    for j = 1:max(size(out))
        y = [y; out{j, i}.mean_time];
        ci = [ci; out{j, i}.CI_time];
    end
    h = plot(vels, y, symbols{i}, vels, ci(:, 1), colors{i}, vels, ci(:, 2), colors{i});
    legend_subset(i) = h(1);
end
hold off
xlabel('velocity [Km/h]');
ylabel('[s]');
legend(legend_subset, 'Custom', 'Random 1', 'Random 2');
%%

%% leftover mem at bs
figure;
title('Average leftover data at BSs')
hold on
grid on
for i = 1:min(size(out))
    y = [];
    for j = 1:max(size(out))
        y = [y; out{j, i}.mean_mem];
    end
    plot(vels, y / 8e6, symbols{i});
end
hold off
xlabel('velocity [Km/h]');
ylabel('[MBytes]');
legend('Custom', 'Random 1', 'Random 2');
%%

% figure;
% hold on
% x_max = 0;
% for i = 1:min(size(out))
%     y = out{4, i}.ue_buff; % 4 is the idx for 100 km/h
%     x = (1:max(size(y)))' * 0.1; %assuming dt=0.1 s
%     if x(end) > x_max
%         x_max = x(end);
%     end
%     plot(x, y / 8e9, '-');  
% end
% plot([1, x_max], ones(2,1) * double(out{4, 1}.ue_max_buff) / 8e9);
% hold off
% xlabel('time [s]');
% ylabel('[MBytes]');
% legend('Random f', 'Custom', 'Random', 'buffer limit');
% title(strcat('UE buffer (with consumption rate of 0.068 Gbps )'));


%% QoS, WARNING: THIS IS ONLY A SKETCH IDEA

%QoS: evaluate average UE buffer usage (ideal should be around 50-75% ???)
out_buff = cell(1,1);
for v = vels
    vel = num2str(v);
    if v < 100
        vel = strcat('0', vel);
    end
    name = strcat('RESULTS//savings_v', vel, '*');
    all_files = dir(name);
    for func_idx = 1:max(size(all_files))      
        saves = load(strcat('RESULTS//', all_files(func_idx).name));
        saves = saves.savings; %due to load() strangeness
        ue_buff = [];
        for i = 1:max(size(saves)) % # of monte carlo iterations per file
            s = saves{i};
            ue_buff = [ue_buff; mean(s.ue_buffer(15:end-15))];%some UE's positions are degenerate due to simulation starting/ending              
        end
        
        idx = (v-70)/10 + 1;
        out_buff{idx, func_idx} = struct;
        out_buff{idx, func_idx}.buff = ue_buff;
    end
end

figure;
title('UE''s buffer average load')
hold on
grid on
max_buff = double(s.ue_max_buffer);
legend_subset = [];
for i = 1:min(size(out))
    y = [];
    ci = [];
    for j = 1:max(size(out))
        y = [y; mean(out_buff{j, i}.buff)/ max_buff];        
        ci = [ci; ConfIntervals(out_buff{j, i}.buff / max_buff)];
    end
    h = plot(vels, y, symbols{i}, vels, ci(:, 1), colors{i}, vels, ci(:, 2), colors{i});
    legend_subset(i) = h(1);
end
plot([vels(1), vels(end)], ones(2,1) * 0.5, '--');
plot([vels(1), vels(end)], ones(2,1))
plot(vels(1), 1.1); % just to show better the max load line
hold off
% grid off
xlabel('velocity [km/h]');
ylabel('Buffer load [%]');
legend(legend_subset, 'Custom', 'Random 1', 'Random 2', 'Ideal load', 'Max load');



figure;
title('QoS averaged for all velocities');
legend_subset = [];
for func = 1:min(size(out))    
    hold on;
    grid on;
    y = [];
    ci = [];
    for vel = 1:max(size(out))
        lost = out{vel,func}.lost_verbose ./ (out{vel,func}.ue_rates * 1e9);
        wait = out{vel,func}.time_verbose;
%         den = out{vel, func}.chunks_verbose/ (0.068 * 1e9);
%         tmp = (wait + lost) ./ (den + wait);
        den = out{vel,func}.sim_duration; %tot simulation duration
        tmp = (wait + lost) ./ den;
        y = [y; mean(tmp)];       
        ci = [ci; ConfIntervals(tmp)];
    end
    y = max(y, 0);
    h = plot(vels, 1 - y, symbols{func}, vels, 1 - ci(:, 1), colors{func}, vels, 1 - ci(:, 2), colors{func});  
    legend_subset(func) = h(1);
    hold off;
    
    if func == 1
        QoS = struct();
        QoS.qos = 1 - y;
        QoS.vels = vels;
        QoS.CI = ci;
        save('RESULTS//QoS.mat', 'QoS');
    end
end
xlabel('Velocity [km/h]');
ylabel('QoS [%]');
legend(legend_subset, 'Custom', 'Random 1', 'Random 2', 'Location', 'southeast');
%%