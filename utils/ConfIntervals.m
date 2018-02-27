function [CI] = ConfIntervals(x, varargin)
    a = 0.05;
    if nargin > 1
        a = varargin{1};
    end                      % Create Data
    SEM = std(x)/sqrt(length(x));               % Standard Error
    ts = tinv([a/2  1-a/2], length(x)-1);      % T-Score
    CI = mean(x) + ts*SEM;

end