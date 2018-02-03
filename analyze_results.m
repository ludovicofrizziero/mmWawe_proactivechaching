% clear all;
% close all;
% clc;
% 
% all_files = dir('RESULTS//savings*');
% load(strcat('RESULTS//', all_files(3).name));
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

out = cell(1,1);
vels = 70:10:130;
for v = vels
    name = strcat('RESULTS//savings_v', num2str(v), '*');
    all_files = dir(name);
    for func_idx = 1:max(size(all_files))      
        saves = load(strcat('RESULTS//', all_files(func_idx).name));
        saves = saves.savings; %due to load() strangeness
        ue_d_lost = [];
        ue_w_time = [];
        bs_mem_left = [];        
        tot_chunks = [];
        for i = 1:max(size(saves)) % # of monte carlo iterations per file
            s = saves{i};
            ue_d_lost = [ue_d_lost; sum(s.ue_lost_data(1:end-5))];%the last UE's positions are degenerate due to simulation ending  
            ue_w_time = [ue_w_time; sum(s.ue_waiting_time(1:end-5))];
            bs_mem_left = [bs_mem_left; sum(s.BSs_mem_state(1:end-15))]; 
            tot_chunks = [tot_chunks; sum(s.chunks)];
        end
        
        idx = (v-70)/10 + 1;
        out{idx, func_idx} = struct;
        out{idx, func_idx}.mean_lost = mean(ue_d_lost);
        out{idx, func_idx}.mean_time = mean(ue_w_time);
        out{idx, func_idx}.mean_mem = mean(bs_mem_left);
        out{idx, func_idx}.time_verbose = ue_w_time;
        out{idx, func_idx}.lost_verbose = ue_d_lost + bs_mem_left;
        out{idx, func_idx}.chunks_verbose = tot_chunks;
        %out{idx, func_idx}.ue_buff = saves{1}.ue_buffer;
        %out{idx, func_idx}.ue_max_buff = saves{1}.ue_max_buffer;
    end
end

%% lost memory
figure;
title('average lost data');
hold on
for i = 1:min(size(out))
    y = [];
    for j = 1:max(size(out))
        y = [y; out{j, i}.mean_lost];
    end
    plot(vels, y / 8e9, '-o');
end
hold off
xlabel('velocity [Km/h]');
ylabel('[GBytes]');
legend('Random f', 'Custom', 'Random');
%%

%% wait time
figure;
title('average cumulative wait time')
hold on
for i = 1:min(size(out))
    y = [];
    for j = 1:max(size(out))
        y = [y; out{j, i}.mean_time];
    end
    plot(vels, y, '-^');
end
hold off
xlabel('velocity [Km/h]');
ylabel('[s]');
legend('Random f', 'Custom', 'Random');
%%

%% leftover mem at bs
figure;
title('average leftover data at bs')
hold on
for i = 1:min(size(out))
    y = [];
    for j = 1:max(size(out))
        y = [y; out{j, i}.mean_mem];
    end
    plot(vels, y / 8e9, '-*');
end
hold off
xlabel('velocity [Km/h]');
ylabel('[GBytes]');
legend('Random f', 'Custom', 'Random');
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
    name = strcat('RESULTS//savings_v', num2str(v), '*');
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
        out_buff{idx, func_idx}.mean = mean(ue_buff);
    end
end

figure;
title('UE buffer average load (with consumption rate 0.068 Gbit/s)')
hold on
grid on
max_buff = double(s.ue_max_buffer);
for i = 1:min(size(out))
    y = [];
    for j = 1:max(size(out))
        y = [y; out_buff{j, i}.mean];
    end
    plot(vels, y / max_buff, '-*');
end
plot([vels(1), vels(end)], ones(2,1) * 0.5, '--');
plot([vels(1), vels(end)], ones(2,1))
plot(vels(1), 1.1); % just to show better the max load line
hold off
grid off
xlabel('velocity [Km/h]');
ylabel('Buffer load [%]');
legend('Random f', 'Custom', 'Random', 'Ideal load', 'Max load');



figure;
title('QoS averaged for all velocities ');
a = 0:0.01:1;
for func = 1:3    
    hold on;
    y = [];
    for vel = 1:max(size(out))
        lost = out{vel,func}.lost_verbose / (0.068 * 1e9);
        wait = out{vel,func}.time_verbose;
        den = out{vel, func}.chunks_verbose/ (0.068 * 1e9);
        tmp = (wait + lost) ./ den;
        y = [y; mean(tmp)];
    end
    plot(vels, 1 - y);    
    hold off;
end
xlabel('Velocity [Km/h]');
ylabel('QoS');
legend('Random f', 'Custom', 'Random');
%%

%% QoS
% for func = 1:3
%     figure;
%     title(strcat('QoS for func ', num2str(func)));
%     a = 0.01:0.01:1;
%     hold on;
%     for vel = 1:max(size(out))
%         lost = out{vel,func}.mean_lost / 0.068 + out{vel, func}.mean_mem / 0.068;
%         wait = out{vel,func}.mean_time;
% 
%         y = (a * 1 / wait + (1-a)* 1 / lost);
%         plot(a, y);
%     end
%     legend('70 Km/h', '80 Km/h', '90 Km/h', '100 Km/h', '110 Km/h', '120 Km/h', '130 Km/h');
%     hold off;
% end
%%