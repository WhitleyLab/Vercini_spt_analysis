% This function fits a sum of three independent Gaussian functions to
% single-particle tracking speeds. It then calculates the expected values
% of each distribution, with errors from bootstrapping.

function [m, res, CI, param] = fit_2gauss(data)

histx=0:2:500;
[bars, edges] = histcounts(data,histx,'Normalization','Probability');
edges = edges(1:end-1) + (edges(2)-edges(1))/2; % shift bin centers by half a bin width

gauss = @(a,x) a(1)*exp(-(x-a(2)).^2 / (2*a(3))); % Gaussian function
gauss2 = @(a,x) gauss([a(1) 0 a(2)],x) + gauss([a(3) a(4) a(5)],x); % Sum of 2 Gaussians, with first fixed at 0
lb = [0 0 0 0 0];
ub = [1 Inf 1 Inf Inf];
initval = [0.5 20 0.002 180 8000];

options = optimoptions('lsqcurvefit','Display','off');
[m, ~, res, ~, ~, ~, J] = lsqcurvefit(gauss2, initval, edges, bars, lb, ub, options);
CI = nlparci(m, res, 'Jacobian', J);

param.hist_xlim = [histx(1) histx(end)];
param.hist_bin_fit = histx(2)-histx(1);
param.fit_fun = gauss2;
param.lower_bounds = lb;
param.upper_bounds = ub;
param.initial_guesses = initval;
param.fit_options = options;

end