%% Function to compute spacial signature for channel matrix H
%
% Frizziero
% 12 / 27 /2017
%
% Giordani - Rebato
% 12/06/2015
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [spatial_sig] = compute_spatial_signature_ssf (K, L, theta, phi, ant_pos, horiz_angle, subpath_angle, lambda)

%% parameters
rms_horiz = 10.2*2*pi/360;                                              % BS cluster rms angular spread for the horizzontal component

%% computation
s = sum(L);
spatial_sig = zeros(max(size(ant_pos)), s);

index = 1;
for k = 1 : K
    if k == 1
        horiz_angle(k) = theta;
%         vertical_angle = phi;
%     else
%         vertical_angle = phi;
    end
    vertical_angle = phi;
    
    for l = 1 : L(k)
        %spatial_sig (index,:) = generate_single_path (horiz_angle(k) + ((-1)^(l-1))*subpath_angle(k,l), vertical_angle, n); % row 455 mmwave-channel-matrix.cc
        spatial_sig (:, index) = steering_vector(horiz_angle(k), vertical_angle, lambda, ant_pos); 
        index = index+1;
    end
    
end

end