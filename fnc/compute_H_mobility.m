%% Test to compute the channel matrix H following ns-3 model
%
% Giordani - Rebato
%
% 12/06/2015
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [H H_params] = compute_H_mobility(f_c,theta_rx,theta_tx,n_rx,n_tx)
%% Parameters
cdf_of_cluster_num = [0.48,0.761,0.927];
xi = 4;
r_tau = 2.8;
phi_rx = 0;                                                                % (Note 3)
phi_tx = 0;                                                                % (Note 3)
BW = 1*10^9;                                                               % BW
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
L = zeros(1 , K);                                                          % number of subpaths per cluster
for k = 1 : K
    L(k) = 1 + round(9 * rand(1,1));
end



%% Computation of power fraction for each cluster

for k = 1 : K
    U_k = rand(1);
    Z_k = xi*randn(1,1);
    cluster_power_fraction(k) = (U_k^(r_tau-1)*10^(0.1*Z_k))/L(k);
end

index = 1;
for k = 1:K
    for l = 1: L(k)
        power_fraction(index) = cluster_power_fraction(k)*10^(0.6*rand(1));
        index = index+1;
    end
end

power_fraction = sort(power_fraction, 'descend');


%% Normalization cluster power fraction

power_sum = 0;
s = sum(L);

for i = 1 : s
    power_sum = power_sum + power_fraction(i);
end

for i = 1 : s
    power_fraction(i) = power_fraction(i) / power_sum;
end


%% Power delay

cluster_dealy = truncated_exprnd(83,300,K);
cluster_dealy = sort(cluster_dealy);

for k = 1:K
    if k == 1
        cluster_dealy(k) = 0;
    else
        cluster_dealy(k) = (cluster_dealy(k-1) + cluster_dealy(k) - cluster_dealy(1) + 25)*(1e-9);
    end
end



%% Subpath dealy
index = 1;
for k = 1:K
    for l=1:L(k)
        subpath_delay(index) = cluster_dealy (k) + (1/BW*(l-1))^(1+0.43*rand(1));
        index = index +1;
    end
end

subpath_delay = sort(subpath_delay);





%% Computation of spatial signatures (matrix, one row for each subpath in each cluster)

[spatial_matrix_tx,horiz_angle_tx,subpath_angle_tx] = compute_spatial_signature (K, L, theta_tx, phi_tx, n_tx);
[spatial_matrix_rx, horiz_angle_rx,subpath_angle_rx] = compute_spatial_signature (K, L, theta_rx, phi_rx, n_rx);

%% Doppler shift
%phi_l = mod(reshape(subpath_angle_rx(find(subpath_angle_rx~=0)),[1,s]) - motion_direction , 2*pi);
%phi_l = zeros(1,s);

%% Computation of small scale fading gain

for i = 1:s
    delay_scf = exp(2 * pi * 1i * f_c * subpath_delay(i));
    doppler_scf =1;
    small_scale_fading(i) = sqrt(power_fraction(i))*doppler_scf/delay_scf;
end

%% Computation of H matrix

H = zeros (n_rx,n_tx);

for i = 1:s
    H = H + (small_scale_fading(i) * spatial_matrix_rx(i,:).' * conj(spatial_matrix_tx(i,:)));
end

H_params = struct('K',K);
H_params = setfield(H_params,'cluster_power_fraction',cluster_power_fraction);
H_params = setfield(H_params,'cluster_dealy',cluster_dealy);
H_params = setfield(H_params,'subpath_delay',subpath_delay);
H_params = setfield(H_params,'spatial_matrix_tx',spatial_matrix_tx);
H_params = setfield(H_params,'spatial_matrix_rx',spatial_matrix_rx);
H_params = setfield(H_params,'power_fraction',power_fraction);
H_params = setfield(H_params,'subpath_angle_tx',subpath_angle_tx);
H_params = setfield(H_params,'subpath_angle_rx',subpath_angle_rx);
H_params = setfield(H_params,'horiz_angle_tx',horiz_angle_tx);
H_params = setfield(H_params,'horiz_angle_rx',horiz_angle_rx);
H_params = setfield(H_params,'L',L);
end