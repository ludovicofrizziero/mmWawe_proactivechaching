%% Test to compute the channel matrix H following ns-3 model
%
% Giordani - Rebato
%
% 12/06/2015
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [H, H_params] = compute_H_mobility(f_c, angles_rx, angles_tx, ant_pos_rx, ant_pos_tx, BW)

%% Parameters
cdf_of_cluster_num = [0.48,0.761,0.927];
xi = 4;
r_tau = 2.8;
phi_rx = angles_rx(2);                                                                % (Note 3)
phi_tx = angles_tx(2);                                                                % (Note 3)
theta_rx = angles_rx(1);
theta_tx = angles_tx(1);
lambda = physconst('lightspeed') / f_c;
%BW = 1*10^9;                                                               % BW
%f_c = 28*10^9;                                                             % central frequency
%speed = 0;                                                                 % Do not consider velocity in this model

%% Computation of clusters

p_ref = rand(1,1);
K = 0;                                                                     % number of clusters
if (p_ref < cdf_of_cluster_num(1))
    K = 1;
elseif (p_ref < cdf_of_cluster_num(2))
    K = 2;
elseif (p_ref < cdf_of_cluster_num(3))
    K = 3;
else
    K = 4;
end



%% Computation of subpaths
% L = zeros(1 , K);                                                          % number of subpaths per cluster
% for k = 1 : K
%     L(k) = 1 + round(9 * rand(1,1));
% end
L = 1 + round(9 * rand(1, K));



%% Computation of power fraction for each cluster

% for k = 1 : K
%     U_k = rand(1);
%     Z_k = xi*randn(1,1);
%     cluster_power_fraction(k) = (U_k^(r_tau-1)*10^(0.1*Z_k))/L(k);
% end
U_k = (rand(1, K)).^(r_tau-1);
Z_k = 10.^(0.1 * (xi * randn(1, K)));
cluster_power_fraction = (U_k .* Z_k) ./ L;

s = sum(L);
index = 1;
power_fraction = rand(1, K * s);
for k = 1:K
    for l = 1: L(k)
        power_fraction(index) = cluster_power_fraction(k)*10^(0.6*power_fraction(index));
        index = index+1;
    end
end

power_fraction = sort(power_fraction, 'descend');


%% Normalization cluster power fraction

% power_sum = 0;
% 
% for i = 1 : s
%     power_sum = power_sum + power_fraction(i);
% end
% 
% for i = 1 : s
%     power_fraction(i) = power_fraction(i) / power_sum;
% end
power_fraction = power_fraction / sum(power_fraction);



%% Power delay

cluster_dealy = truncated_exprnd(83,300,K);
cluster_dealy = sort(cluster_dealy);

% for k = 1:K
%     if k == 1
%         cluster_dealy(k) = 0;
%     else
%         cluster_dealy(k) = (cluster_dealy(k-1) + cluster_dealy(k) - cluster_dealy(1) + 25)*(1e-9);
%     end
% end
cluster_dealy(1) = 0;
cluster_dealy(2:end) = (cluster_dealy(1:end-1) + cluster_dealy(2:end) + 25) * (1e-9);


%% Subpath dealy
index = 1;
T = 1/BW;
subpath_delay = rand(1, K * s);
for k = 1:K
    for l=1:L(k)
        subpath_delay(index) = cluster_dealy (k) + (T *(l-1))^(1+0.43*subpath_delay(index));
        index = index +1;
    end
end

subpath_delay = sort(subpath_delay);





%% Computation of spatial signatures (matrix, one row for each subpath in each cluster)

[spatial_matrix_tx, horiz_angle_tx, subpath_angle_tx] = compute_spatial_signature (K, L, theta_tx, phi_tx, ant_pos_tx, lambda);
[spatial_matrix_rx, horiz_angle_rx, subpath_angle_rx] = compute_spatial_signature (K, L, theta_rx, phi_rx, ant_pos_rx, lambda);

%% Doppler shift
%phi_l = mod(reshape(subpath_angle_rx(find(subpath_angle_rx~=0)),[1,s]) - motion_direction , 2*pi);
%phi_l = zeros(1,s);

%% Computation of small scale fading gain

% for i = 1:s
%     delay_scf = exp(2 * pi * 1i * f_c * subpath_delay(i));
%     doppler_scf =1;
%     small_scale_fading(i) = sqrt(power_fraction(i))*doppler_scf/delay_scf;
% end
delay_scf = exp((2 * pi * 1i * f_c) * subpath_delay);
small_scale_fading = sqrt(power_fraction) ./ delay_scf;


%% Computation of H matrix

H = zeros (max(size(ant_pos_rx)), max(size(ant_pos_tx)));

for i = 1:s
    H = H + (small_scale_fading(i) * spatial_matrix_rx(:, i) * spatial_matrix_tx(:, i)');
end

H_params = struct('K',K);
H_params.('cluster_power_fraction') = cluster_power_fraction;
H_params.('cluster_dealy') = cluster_dealy;
H_params.('subpath_delay') = subpath_delay;
H_params.('spatial_matrix_tx') = spatial_matrix_tx;
H_params.('spatial_matrix_rx') = spatial_matrix_rx;
H_params.('power_fraction') = power_fraction;
H_params.('subpath_angle_tx') = subpath_angle_tx;
H_params.('subpath_angle_rx') = subpath_angle_rx;
H_params.('horiz_angle_tx') = horiz_angle_tx;
H_params.('horiz_angle_rx') = horiz_angle_rx;
H_params.('L') = L;

end