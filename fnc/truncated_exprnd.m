%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Assume that lambda is the exponential parameter, and T the truncation point.
% Generate TRUNCATED EXPONENTIAL RV
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [truncated_variables] = truncated_exprnd (mu, bound , n)

  R = rand(n,1)*(1-exp(-bound/mu));
  truncated_variables = -log(1-R)*mu;
end