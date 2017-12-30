%% Test to compute the channel matrix H following ns-3 model
%
% In this case, I just have to recompute the small-scale fading (SSF),
% since it is assumed that the large scale fading (angular spread of
% cluster) hasn't changed in this time interval.
%
% ASSUMPTION: in this version of the code, just change the Doppler
% parameter, all the other parameters are kept from previous computations
% of the H matrix.
%
% Giordani - Rebato
%
% 12/06/2015
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [H] = compute_H_ssf(f_c, angles_rx, angles_tx, ant_pos_rx, ant_pos_tx, H_params, BW)


%% Parameters
cdf_of_cluster_num = [0.48,0.761,0.927];
xi = 4;
r_tau = 2.8;
phi_rx = angles_rx(2);                                                                % (Note 3)
phi_tx = angles_tx(2);                                                                % (Note 3)
theta_rx = angles_rx(1);
theta_tx = angles_tx(1);
%BW = 1*10^9;                                                               % BW
%f_c = 28*10^9;                                                             % central frequency
%speed = 0;                                                                 % Do not consider velocity in this model
lambda = physconst('lightspeed') / f_c;

%% Computation of clusters (is used in the ssc

K = H_params.K; %take K from previous case
 
L = H_params.L; %take L from previous case

power_fraction = H_params.power_fraction; %take power_fraction from previous case, it should be already normalized.

s = sum(L);


%% Power delay

%cluster_dealy = H_params.cluster_dealy; %take cluster delays from previous case

%% Subpath dealy
% index = 1;
% for k = 1:K
%     for l=1:L(k)
%         subpath_delay(index) = cluster_dealy (k) + (1/BW*(l-1))^(1+0.43*rand(1));
%         index = index +1;
%     end
% end
% 
% subpath_delay = sort(subpath_delay);
subpath_delay = H_params.subpath_delay;


%% Doppler shift
%phi_l = zeros(1,s);

%% Computation of spatial signatures (matrix, one col for each subpath in each cluster)

%take spatial signatures from previous case
%spatial_matrix_tx = H_params.spatial_matrix_tx;
%spatial_matrix_rx = H_params.spatial_matrix_rx;
spatial_matrix_tx = compute_spatial_signature_ssf (H_params.K, H_params.L, theta_tx, phi_tx, ant_pos_tx,H_params.horiz_angle_tx,H_params.subpath_angle_tx, lambda);
spatial_matrix_rx = compute_spatial_signature_ssf (H_params.K, H_params.L, theta_rx, phi_rx, ant_pos_rx,H_params.horiz_angle_rx,H_params.subpath_angle_rx, lambda);

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
        H = H + ( small_scale_fading(i) * spatial_matrix_rx(:, i) * spatial_matrix_tx(:, i)' );
end

end