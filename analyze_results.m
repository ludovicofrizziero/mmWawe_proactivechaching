clear all;
close all;
clc;

load('RESULTS//savings.mat');

s = savings{1};

figure;
plot(s.servingBS_IDs, '-');
title('serving BS ID');

figure;
hold on
plot(1:max(size(s.BSs_mem_state)), s.chunks, 'o');
plot(1:max(size(s.BSs_mem_state)), s.BSs_mem_state, '*');
title('memory state at BSs');
legend('starting situation', 'final situation');
hold off

figure;
plot(s.rate, '-');
title('rate [bit/s]');

figure;
plot(s.ue_buffer);
title(strcat('UE buffer (with consumption rate of', num2str(s.ue_requested_rate, ' %1.3f'), ' Gbps )'));

figure;
plot(s.ue_waiting_time, 'o');
title('UE waiting times');

figure;
plot(s.ue_lost_data);
title('UE lost data');



