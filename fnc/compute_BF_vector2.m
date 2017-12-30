%%
% Frizziero
% AoA: angle of arrival of the wawe.
% lambda: wawelength.
% ant_pos: array [M x 1] containing the (x,y,z) position of each antenna (that are M in total)
%       expressed in wawelengths.
% inv_cov_mat: inverse of the spatial covariance matrix (M x M).
%%
function [w] = compute_BF_vector2(angles, lambda, ant_pos)
    sigma_noise = 0.01;
    theta = angles(1); %azimut
    phi = angles(2); %elevation
    v = steering_vector(theta, phi, lambda, ant_pos);
    
    R = v * v' + sigma_noise^2 * eye(max(size(v))); %cov matrix, assuming signal power of 1 and no interference
    tmp = R \ v;
    w = tmp / (v' * tmp); %MVDR method [R^-1 * v * (v^H * R^-1 * v)^-1]
    
    %w' * steer_vec %array response. Should be 1+0i
    
    %% possible extension:
    % w = inv(R) * C / (C' * inv(R) * C) * F
    % where f is a gain vector with ones in interesting direction and zeros
    % in interfering direction or something appropriate anyway [dim = (Q x 1)]
    % C is a matrix of all steering vectors [dim = (M x Q)]   
end