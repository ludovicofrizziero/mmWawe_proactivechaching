%% Code to compute the BF vector
%
% Giordani - Rebato
%
% 11/19/2015
%
% Note 1: vertical angle is fixed to zero.
%

function [BF_vector] = compute_BF_vector(n_elements,angle)

%% Parameters
phi = 0;
theta = angle;


%% Computation
index = 1;
BF_vector = zeros(1,n_elements);
for vIndex = 1 : sqrt(n_elements)
    for hIndex = 1 : sqrt(n_elements)
        BF_vector(index) = sqrt(1/n_elements)*exp(1i*pi*(hIndex-1)*sin(theta) + 1i*pi*(vIndex-1)*sin(phi));
        
%         BF_vector(index) = sqrt(1/n_elements)*exp(1i*(hIndex-1)*(theta) + 1i*(vIndex-1)*(phi)); % BEAMSPACED BF VECTORS

        index=index+1;
    end
end

end