%% Function to compute spacial signature for channel matrix H
%
% Giordani - Rebato
%
% 12/06/2015
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [spatial_sig] = compute_spatial_signature_ssf (K, L, theta, phi, n, horiz_angle, subpath_angle)

%% parameters
rms_horiz = 10.2*2*pi/360;                                              % BS cluster rms angular spread for the horizzontal component

%% computation
s = sum(L);
spatial_sig = zeros(s,n);

index = 1;
for k = 1 : K
    if k == 1
        horiz_angle(k) = theta;
        vertical_angle = phi;
    else
        vertical_angle = phi;
    end
    
    for l = 1 : L(k)
        %spatial_sig (index,:) = generate_single_path (horiz_angle(k) + ((-1)^(l-1))*subpath_angle(k,l), vertical_angle, n); % row 455 mmwave-channel-matrix.cc
        spatial_sig (index,:) = generate_single_path (horiz_angle(k), vertical_angle, n); % row 455 mmwave-channel-matrix.cc
        index = index+1;
    end
    
end

end