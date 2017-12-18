function [ PL_LTE ] = PL_LTE( R )
% compute PL of UE at distance R wrt the LTE BS
% OBS: R is in km

variable = rand(1); % random variable
P_LOS = min(0.018/R,1)*(1-exp(-R/0.063))+exp(-R/0.063); % Probability that the test UE is in LOS

if (variable <= P_LOS) % I am in LOS
    PL_LTE = 103.4 + 24.2*log10(R);
else % NLOS
    PL_LTE = 131.1 + 42.8*log10(R);
end

end

