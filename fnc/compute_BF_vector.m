%% Code to compute the BF vector
%
% Giordani - Rebato
%
% 11/19/2015
%
% Note 1: vertical angle is fixed to zero, if not specified as optional
% parameter
%

function [BF_vector] = compute_BF_vector(n_elements, hor_angle, varargin)

%% Parameters
phi = 0; %elevation of rx???
if nargin > 2
    phi = varargin{1};
end
theta = hor_angle; %horizontal direction of rx

%% precomputed values for speed
BF_vector = zeros(1,n_elements);
n_sqr = sqrt(n_elements);
a = sqrt(1/n_elements);
i_pi_stheta = 1i * pi * sin(theta);
i_pi_sphi = 1i * pi * sin(phi);

%% Computation
index = 1;
for vIndex = 1 : n_sqr
    tmp = (vIndex - 1) * i_pi_sphi;
    for hIndex = 1 : n_sqr
        BF_vector(index) = a * exp(i_pi_stheta * (hIndex-1) + tmp);
      
%         BF_vector(index) = sqrt(1/n_elements)*exp(1i*(hIndex-1)*sin(theta) + 1i*(vIndex-1)*sin(phi)); % BEAMSPACED BF VECTORS

        index=index+1;
    end
end

end