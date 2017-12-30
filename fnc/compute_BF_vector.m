%% Code to compute the BF vector
%
% Giordani - Rebato - Frizziero
%
% 12/27/2017
%
%
%

function [BF_vector] = compute_BF_vector(n_elements, hor_angle, vert_angle)

%% Parameters
phi = vert_angle;
theta = hor_angle; 


%% precomputed values for speed
BF_vector = zeros(1,n_elements);
n_sqr = sqrt(n_elements);
a = sqrt(1/n_elements);
i_pi_stheta = 1i * pi * sin(theta) * d;
i_pi_sphi = 1i * pi * sin(phi) * d;

%% Computation
index = 1;
for vIndex = 1 : n_sqr
    tmp = (vIndex - 1) * i_pi_sphi;
    for hIndex = 1 : n_sqr
          BF_vector(index) = a * exp(i_pi_stheta * (hIndex-1) + tmp);

%         BF_vector(index) = a * exp(1i*(hIndex-1)*sin(theta) + 1i*(vIndex-1)*sin(phi)); % BEAMSPACED BF VECTORS

        index=index+1;
    end
end   

end