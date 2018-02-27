%% 
% N : number of samples
% ant_array : number of antennas per BS (optional, default is 1)
% users : number of users the antenna is serving (optional, default is 1)
% mod : modulation scheme to use, either 'DEF' or ... [not implemented yet] (optional, default is 'DEF' with unit power signals) 
% power : noise power (optional, default is 0)
% var : noise variance (optional, default is 1)
%% 
function [samples] = generate_signal_samples(N, varargin)
    
    users = 1;
    ant_array = 1;
    mod = 'DEF';
    n_var = 1;
    s_power = 0;
    n_power = 0;
    
    switch nargin
        case 2
            ant_array = varargin{2};
        case 3
            ant_array = varargin{2};
            users = varargin{3};
        case 4
            ant_array = varargin{2};
            users = varargin{3};
            mod = varargin{4};
        case 5
            ant_array = varargin{2};
            users = varargin{3};
            mod = varargin{4};
            n_power = varargin{5};
        case 6
            ant_array = varargin{2};
            users = varargin{3};
            mod = varargin{4};
            n_power = varargin{5};
            n_var = varargin{6};
    end 
    
    switch mod
        case 'DEF'
            s_power = 1;
        case 'PAM'
            s_power = 1; %not implemented yet
        case 'QAM'
            s_power = 1; %not implemented yet
    end
    
    samples = cell(N);
    for i = 1:N
        s = s_power * ones(ant_array, users); %signal
        interf = 0 * ones(ant_array, users); %interference
        noise = n_power + sqrt(n_var/2) * (randn(ant_array, users) + 1i*randn(ant_array, users)); %complex symmetrical gaussian noise of (0, var)
        samples{i} = s + interf + noise;
    end
    
end