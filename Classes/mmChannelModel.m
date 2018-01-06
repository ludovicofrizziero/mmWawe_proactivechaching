%% 
% Frizziero 1/2/2018
%%

classdef mmChannelModel < handle
    properties (Access = private)
        H
        H_params
        BW
        f
        xi
        r_tau
                
        lambda
    end
    
    methods
        function self =  mmChannelModel(BW, f, xi, r_tau)
            self.BW = BW;
            self.f = f;           
            self.lambda = physconst('lightspeed') / f;
            self.xi = xi; %4;
            self.r_tau = r_tau; %2.8;
        end
        
        function update_channel_state(self, angles_rx, angles_tx, ant_pos_rx, ant_pos_tx, update_all_params)
            if update_all_params
                self.compute_H(angles_rx, angles_tx, ant_pos_rx, ant_pos_tx);
            else
                self.compute_H_ssf(angles_rx, angles_tx, ant_pos_rx, ant_pos_tx);
            end          
        end       
        
        function out = mtimes(a, b) %overload matrix multiplication operator
            try
                out = a.H * b; %a is a mmChannelModel object, b another matrix
            catch E
                out = a * b.H; %b is a mmChannelModel object, a another matrix
            end
        end
    end
    
    methods (Access = private)
        function compute_H(self, angles_rx, angles_tx, ant_pos_rx, ant_pos_tx)
            phi_rx = angles_rx(2);                                                                
            phi_tx = angles_tx(2);                                                                
            theta_rx = angles_rx(1);
            theta_tx = angles_tx(1);
            
            %% Computation of clusters
            p_ref = rand(1,1);
            cdf_of_cluster_num = [0.48, 0.761, 0.927];
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
            L = 1 + round(9 * rand(1, K));
            
            %% Computation of power fraction for each cluster
            U_k = (rand(1, K)).^(self.r_tau-1);
            Z_k = 10.^(0.1 * (self.xi * randn(1, K)));
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
            power_fraction = power_fraction / sum(power_fraction);
            
            %% Power delay
            cluster_dealy = truncated_exprnd(83,300,K);
            cluster_dealy = sort(cluster_dealy);
            cluster_dealy(1) = 0;
            cluster_dealy(2:end) = (cluster_dealy(1:end-1) + cluster_dealy(2:end) + 25) * (1e-9);
            
            %% Subpath dealy
            index = 1;
            T = 1 / self.BW;
            subpath_delay = rand(1, K * s);
            for k = 1:K
                for l=1:L(k)
                    subpath_delay(index) = cluster_dealy (k) + (T *(l-1))^(1+0.43*subpath_delay(index));
                    index = index +1;
                end
            end

            subpath_delay = sort(subpath_delay);
            
            %% spatial signature
            [spatial_matrix_tx, horiz_angle_tx, subpath_angle_tx] = compute_spatial_signature (K, L, theta_tx, phi_tx, ...
                                                                                                    ant_pos_tx, self.lambda);
            [spatial_matrix_rx, horiz_angle_rx, subpath_angle_rx] = compute_spatial_signature (K, L, theta_rx, phi_rx, ...
                                                                                                    ant_pos_rx, self.lambda);

            %% Computation of small scale fading gain

            % for i = 1:s
            %     delay_scf = exp(2 * pi * 1i * f_c * subpath_delay(i));
            %     doppler_scf =1;
            %     small_scale_fading(i) = sqrt(power_fraction(i))*doppler_scf/delay_scf;
            % end
            delay_scf = exp((2 * pi * 1i * self.f) * subpath_delay);
            small_scale_fading = sqrt(power_fraction) ./ delay_scf;
            
            %% Computation of H matrix
%             H_ = zeros(max(size(ant_pos_rx)), max(size(ant_pos_tx)));
% 
%             for i = 1:s
%                 H_ = H_ + (small_scale_fading(i) * spatial_matrix_rx(:, i) * spatial_matrix_tx(:, i)');
%             end            
%             self.H = H_; %for speed

            D = diag(small_scale_fading(1:s));
            self.H = spatial_matrix_rx * D * spatial_matrix_tx';

            H_params_ = struct('K', K);
            H_params_.('cluster_power_fraction') = cluster_power_fraction;
            H_params_.('cluster_dealy') = cluster_dealy;
            H_params_.('subpath_delay') = subpath_delay;
            H_params_.('spatial_matrix_tx') = spatial_matrix_tx;
            H_params_.('spatial_matrix_rx') = spatial_matrix_rx;
            H_params_.('power_fraction') = power_fraction;
            H_params_.('subpath_angle_tx') = subpath_angle_tx;
            H_params_.('subpath_angle_rx') = subpath_angle_rx;
            H_params_.('horiz_angle_tx') = horiz_angle_tx;
            H_params_.('horiz_angle_rx') = horiz_angle_rx;
            H_params_.('L') = L;
            self.H_params = H_params_; %for speed
            
        end
        
        function compute_H_ssf(self, angles_rx, angles_tx, ant_pos_rx, ant_pos_tx) %small scale fading computation only
            %% Parameters
            phi_rx = angles_rx(2);                                                                
            phi_tx = angles_tx(2);                                                                
            theta_rx = angles_rx(1);
            theta_tx = angles_tx(1);

            %% Computation of clusters (is used in the ssc

            K = self.H_params.K; %take K from previous case

            L = self.H_params.L; %take L from previous case

            power_fraction = self.H_params.power_fraction; %take power_fraction from previous case, it should be already normalized.

            s = sum(L);
           
            %% Subpath dealy            
            subpath_delay = self.H_params.subpath_delay;

            spatial_matrix_tx = compute_spatial_signature_ssf (self.H_params.K, self.H_params.L, theta_tx, phi_tx, ...
                                            ant_pos_tx, self.H_params.horiz_angle_tx, self.H_params.subpath_angle_tx, self.lambda);
            spatial_matrix_rx = compute_spatial_signature_ssf (self.H_params.K, self.H_params.L, theta_rx, phi_rx, ...
                                            ant_pos_rx, self.H_params.horiz_angle_rx, self.H_params.subpath_angle_rx, self.lambda);

            %% Computation of small scale fading gain
            delay_scf = exp((2 * pi * 1i * self.f) * subpath_delay);
            small_scale_fading = sqrt(power_fraction) ./ delay_scf;

            %% Computation of H matrix
%             self.H = zeros (max(size(ant_pos_rx)), max(size(ant_pos_tx)));
% 
%             for i = 1:s
%                     self.H = self.H + ( small_scale_fading(i) * spatial_matrix_rx(:, i) * spatial_matrix_tx(:, i)' );
%             end

            D = diag(small_scale_fading(1:s));            
            self.H = spatial_matrix_rx * D * spatial_matrix_tx';
        end
    end
end