function [pproc] = poisson2d_fixed(lambda)
% the number of points is Poisson(lambda)-distributed



% conditioned that the number of points is N,
% the points are uniformly distributed
pproc = rand(lambda, 2);
