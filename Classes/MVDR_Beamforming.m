%% 
% Frizziero 1/2/2018
%%

classdef MVDR_Beamforming < handle
    properties (Access = private)
        ant_pos %antenna array fixed structure
        white_noise
        lambda        
        w %beamforming vector
    end
    
    methods
        function self = MVDR_Beamforming(white_noise, lambda, ant_pos)
            self.ant_pos = ant_pos;
            self.white_noise = white_noise;      
            self.lambda = lambda;            
        end
        
        %%
        % angles: [azimut, elevation] of the desired signal arrival (or departure) direction
        % varargin: optional array N x 3 interference arrival directions and relative interference power 
        %                                                               ([interf_power_stdvar, azimut, elevation])
        %%
        function update_state(self, angles, varargin)
            v = steering_vector(angles(1), angles(2), self.lambda, self.ant_pos); 
            R = v * v' + self.white_noise^2 * eye(max(size(v))); %assume unit power for desired signal, then add noise
            
            if nargin == 3 % nargin > 3 --> we have an error!!!
                dirs = varargin{3};      
                for i = 1:max(size(dirs))
                    tmp = steering_vector(dirs(i, 2), dirs(i, 3), self.lambda, self.ant_pos);
                    R = R + dirs(i, 1)^2 + tmp * tmp'; %interference contribution
                end
            end
            
            tmp = R \ v;
            w_tmp = tmp / (v' * tmp); %MVDR method [R^-1 * v * (v^H * R^-1 * v)^-1]   
            self.w = w_tmp / norm(w_tmp);
            
            %w' * v  %array response. Should be 1+0i, check for debug if needed

            %% possible extension:
            % w = inv(R) * C / (C' * inv(R) * C) * F
            % where f is a gain vector with ones in interesting direction and zeros
            % in interfering direction or something appropriate anyway [dim = (Q x 1)]
            % C is a matrix of all steering vectors [dim = (M x Q)]
        end
        
        function out = mtimes(a, b) %overload matrix multiplication            
            try
                out = a * b.w; %a is a matrix N x length(w) [N >= 1], and b is a MVDR object
            catch E
                out = a.w * b; % b is a matrix length(w) x N [N >= 1], and a is a MVDR object
            end %if neither of the two cases where correct, then another exception is trown since we have an incorrect matrix/vector multiplication
        end
        
        function out = ctranspose(a) %overload complex conjugate transpose
            out = a.w';
        end
        
        function varargout = size(this, varargin) %overload size operator
            [varargout{1:nargout}] = builtin('size', this.w, varargin{:});
        end
    end
       
%     methods (Static)
%         function v = steeringVector(theta, phi, lambda, ant_pos)
%             k = (2*pi/lambda) * [sin(phi)*cos(theta); sin(phi)*sin(theta); cos(phi)]; %wawevector
%             v = exp(-1i* ant_pos * k); %should be a length(arr_pos) x 1 vector
%         end
%     end
end