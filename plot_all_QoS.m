clear all;
close all;

name ='RESULTS//QoS*';
all_files = dir(name);

symbols = {'b-*', 'r-^', 'k-o', 'g-x'}; %for functions
colors = {'b:', 'r:', 'k:', 'g:'}; %for Conf. Interv.

figure;
hold on
grid on
% WARNING: IT IS ASSUMED THAT FILES ARE ORDERED ALPHABETICALLY WITH BS DENSITY
for i = 1:length(all_files)
    file = load(strcat('RESULTS//', all_files(i).name));
    QoS = file.QoS;
    h = plot(QoS.vels, QoS.qos, symbols{i}, QoS.vels, 1 - QoS.CI(:, 1), colors{i}, QoS.vels, 1 - QoS.CI(:, 2), colors{i});
    legend_subset(i) = h(1);
end
title('Quality of Service for different BS densities');
xlabel('velocity [Km/h]');
ylabel('QoS [%]');
legend(legend_subset, '10 BS per km', '14 BS per km', '20 BS per km', 'Location', 'southeast');
hold off