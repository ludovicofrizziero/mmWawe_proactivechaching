%Frizziero 27/12/2017

function [v] = steering_vector(theta, phi, lambda, ant_pos)
    k = (2*pi/lambda) * [sin(phi)*cos(theta); sin(phi)*sin(theta); cos(phi)]; %wawevector
    v = exp(-1i* ant_pos * k); %should be a length(arr_pos) x 1 vector
end