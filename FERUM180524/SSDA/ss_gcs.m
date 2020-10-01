function subtempu = ss_gcs(subgermU,rho,rand_generator)

% Generate a random walk from the current germ simu using Gaussian
% conditional sampling
%
% subtempu: generated samples

% Generate realizations of random U-vector
switch rand_generator
    case 0  % Matlab function
        subtempu = randn(size(subgermU));
    otherwise  % twister
        subtempu = inv_norm_cdf(twister(size(subgermU)));
end

subtempu = diag(rho)*subgermU + diag(sqrt(1-rho.^2))*subtempu;  % one move of Markov Chain

      
