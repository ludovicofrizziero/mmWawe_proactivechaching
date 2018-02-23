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
    name = strcat('RESULTS//savings20km_v', vel, '*');
    all_files = dir(name);
    for func_idx = 1:max(size(all_files))      
        saves = load(strcat('RESULTS//', all_files(func_idx).name));
        saves = saves.savings; %due to load() strangeness
        dims = size(saves);
        ue_d_lost = zeros(dims(1), 1);
        ue_w_time = zeros(dims(1), 1);
        bs_mem_left = zeros(dims(1), 1);        
        tot_chunks = zeros(dims(1), 1);
        bs_rates = cell(dims(1), 1);
        ue_rates = zeros(dims(1), 1);
        sim_duration = zeros(dims(1), 1);        
        for itr = 1:dims(1) % # of monte carlo iterations per file
            bs_rates_tmp = [];
            for km = 1:dims(2)
                s = saves{itr, km};
                ue_d_lost(itr) = ue_d_lost(itr) + sum(s.ue_lost_data);%the last UE's positions are degenerate due to simulation ending  
                ue_w_time(itr) = ue_w_time(itr) + sum(s.ue_waiting_time);
                bs_mem_left(itr) = bs_mem_left(itr) + sum(s.BSs_mem_state); 
                tot_chunks(itr) = tot_chunks(itr) + sum(s.chunks);
                ue_rates(itr) = s.ue_requested_rate;
                sim_duration(itr) = sim_duration(itr) + length(s.ue_waiting_time) * 0.1;
                bs_rates_tmp = [bs_rates_tmp; s.rate];
                if sum(s.chunks) == 0
                    tot_chunks(end) = 1e9; %sometimes a solution is not found since all BS have full memory
                end
            end
            bs_rates{itr} = bs_rates_tmp;
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
        out{idx, func_idx}.sim_duration = sim_duration;
        out{idx, func_idx}.ue_rates = ue_rates;
        %out{idx, func_idx}.ue_buff = saves{1}.ue_buffer;
        %out{idx, func_idx}.ue_max_buff = saves{1}.ue_max_buffer;
    end
end

symbols = {'b-*', 'r-^', 'k-o', 'g-x'}; %for functions
colors = {'b:', 'r:', 'k:', 'g:'}; %for Conf. Interv.

%% lost memory
figure;
title('average lost data');
hold on
grid on
for itr = 1:min(size(out))
    y = [];
    for j = 1:max(size(out))
        y = [y; out{j, itr}.mean_lost];
    end
    plot(vels, y / 8e6, symbols{itr});
end
hold off
xlabel('velocity [Km/h]');
ylabel('[MBytes]');
legend('Custom', 'Random 1', 'Random 2');
%%

%% wait time
figure;
title('average cumulative wait time')
hold on
grid on
legend_subset = [];
for itr = 1:min(size(out))
    y = [];
    ci = [];
    for j = 1:max(size(out))
        y = [y; out{j, itr}.mean_time];
        ci = [ci; out{j, itr}.CI_time];
    end
    h = plot(vels, y, symbols{itr}, vels, ci(:, 1), colors{itr}, vels, ci(:, 2), colors{itr});
    legend_subset(itr) = h(1);
end
hold off
xlabel('velocity [Km/h]');
ylabel('[s]');
legend(legend_subset, 'Custom', 'Random 1', 'Random 2');
%%

%% leftover mem at bs
figure;
title('average leftover data at bs')
hold on
grid on
for itr = 1:min(size(out))
    y = [];
    for j = 1:max(size(out))
        y = [y; out{j, itr}.mean_mem];
    end
    plot(vels, y / 8e6, symbols{itr});
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

% %QoS: evaluate average UE buffer usage (ideal should be around 50-75% ???)
% out_buff = cell(1,1);
% for v = vels
%     vel = num2str(v);
%     if v < 100
%         vel = strcat('0', vel);
%     end
%     name = strcat('RESULTS//savings20km_v', vel, '*');
%     all_files = dir(name);
%     for func_idx = 1:max(size(all_files))      
%         saves = load(strcat('RESULTS//', all_files(func_idx).name));
%         saves = saves.savings; %due to load() strangeness
%         ue_buff = [];
%         for itr = 1:max(size(saves)) % # of monte carlo iterations per file
%             s = saves{itr};
%             ue_buff = [ue_buff; mean(s.ue_buffer(15:end-15))];%some UE's positions are degenerate due to simulation starting/ending              
%         end
%         
%         idx = (v-70)/10 + 1;
%         out_buff{idx, func_idx} = struct;
%         out_buff{idx, func_idx}.buff = ue_buff;
%     end
% end
% 
% figure;
% title('UE buffer average load')
% hold on
% grid on
% max_buff = double(s.ue_max_buffer);
% legend_subset = [];
% for itr = 1:min(size(out))
%     y = [];
%     ci = [];
%     for j = 1:max(size(out))
%         y = [y; mean(out_buff{j, itr}.buff)/ max_buff];        
%         ci = [ci; ConfIntervals(out_buff{j, itr}.buff / max_buff)];
%     end
%     h = plot(vels, y, symbols{itr}, vels, ci(:, 1), colors{itr}, vels, ci(:, 2), colors{itr});
%     legend_subset(itr) = h(1);
% end
% plot([vels(1), vels(end)], ones(2,1) * 0.5, '--');
% plot([vels(1), vels(end)], ones(2,1))
% plot(vels(1), 1.1); % just to show better the max load line
% hold off
% % grid off
% xlabel('velocity [Km/h]');
% ylabel('Buffer load [%]');
% legend(legend_subset, 'Custom', 'Random1', 'Random2', 'Ideal load', 'Max load');



figure;
title('QoS for 20km averaged for all velocities');
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
end
xlabel('velocity [Km/h]');
ylabel('QoS [%]');
legend(legend_subset, 'Custom', 'Random 1', 'Random 2', 'Location', 'southeast');
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

figure;
plot(bs_rates{1})