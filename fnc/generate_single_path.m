%% Function to compute single path (single vector) for spatial signature
%
% Giordani - Rebato
%
% 12/07/2015
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



function [spatial_signature] = generate_single_path (hor_angle, vert_angle, n)
index = 1;
% for j = 1 : sqrt(n)
%     for k = 1 : sqrt(n)
%         spatial_signature(index) = exp(1i*pi*sin(hor_angle)*(k-1) + 1i*pi*sin(vert_angle)*(j-1));
%         index = index + 1;
%     end
% end

%precompute for speed, added by Frizziero
spatial_signature = zeros(1, n); 
a = 1i*pi*sin(hor_angle);
b = 1i*pi*sin(vert_angle);
s = sqrt(n);

for j = 1 : s
    tmp = b*(j-1);
    for k = 1 : s
        spatial_signature(index) = exp(a*(k-1) + tmp); 
        index = index + 1;
    end
end

end