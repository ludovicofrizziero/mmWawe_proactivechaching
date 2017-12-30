%% Function to compute spacial signature for channel matrix H
%
% Giordani - Rebato
%
% 12/06/2015
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [spatial_sig, horiz_angle, subpath_angle] = compute_spatial_signature (K, L, theta, phi, ant_pos, lambda)

% %% parameters
% rms_horiz = 10.2*2*pi/360;                                              % BS cluster rms angular spread for the horizzontal component
% 
% %% computation
% s = sum(L);
% spatial_sig = zeros(s,n);
% 
% index = 1;
% for k = 1 : K
%     if k == 1
%         horiz_angle(k) = theta;
%         vertical_angle = phi;
%     else
%         horiz_angle(k) = -pi + (2*pi)*rand(1);
%         vertical_angle = phi;
%     end
%     
%     for l = 1 : L(k)
%         subpath_angle(k,l) = truncated_exprnd(rms_horiz,0.7,1)/2;
%         spatial_sig (index,:) = generate_single_path (horiz_angle(k) + ((-1)^(l-1))*subpath_angle(k,l), vertical_angle, n); % row 455 mmwave-channel-matrix.cc
%         index = index+1;
%     end
%     
% end

%% parameters
rms_horiz = 10.2*2*pi/360;                                              % BS cluster rms angular spread for the horizzontal component

%% computation
s = sum(L);
spatial_sig = zeros(max(size(ant_pos)), s);

horiz_angle = zeros(1, K); %added by Frizziero
subpath_angle = zeros(K, max(L)); %added by Frizziero

index = 1;
for k = 1 : K
    if k == 1
        horiz_angle(k) = theta;
        vertical_angle = phi;
    else
        horiz_angle(k) = -pi + (2*pi)*rand(1);
        vertical_angle = phi;
    end
       
    for l = 1 : L(k)
        subpath_angle(k,l) = truncated_exprnd(rms_horiz, 0.7, 1)/2;
        spatial_sig (:, index) = steering_vector(horiz_angle(k) + ((-1)^(l-1))*subpath_angle(k,l), vertical_angle, lambda, ant_pos); % row 455 mmwave-channel-matrix.cc
        index = index+1;
    end
    
end

end