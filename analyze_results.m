clear all;
close all;
clc;

all_files = dir('RESULTS//savings*');
load(strcat('RESULTS//', all_files(1).name));

s = savings{1};

figure;
plot(s.servingBS_IDs, '-');
title('serving BS ID');

figure;
hold on
plot(s.all_ids, s.chunks, 'o');
plot(s.all_ids, s.BSs_mem_state, '*');
title('memory state at BSs');
legend('starting situation', 'final situation');
hold off

figure;
plot(s.rate, '-');
title('rate [bit/s]');

figure;
hold on
plot(s.ue_buffer);
plot([1,max(size(s.ue_buffer))], ones(2,1) * double(s.ue_max_buffer));
hold off
title(strcat('UE buffer (with consumption rate of', num2str(s.ue_requested_rate, ' %1.3f'), ' Gbps )'));

figure;
plot(s.ue_waiting_time, 'o');
title('UE waiting times');

figure;
plot(s.ue_lost_data);
title('UE lost data');



