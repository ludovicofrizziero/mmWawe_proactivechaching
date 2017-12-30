%% Function to compute spacial signature for channel matrix H
%
% Frizziero
% 12 / 27 /2017
%
% Giordani - Rebato
% 12/06/2015
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [spatial_sig] = compute_spatial_signature_fixed (theta, phi, n)

horiz_angle = theta;
vertical_angle = phi;

spatial_sig  = generate_single_path (horiz_angle, vertical_angle, n); % row 455 mmwave-channel-matrix.cc


end