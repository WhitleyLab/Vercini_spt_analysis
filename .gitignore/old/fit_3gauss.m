% This function fits a sum of three independent Gaussian functions to
% single-particle tracking speeds. It then calculates the expected values
% of each distribution, with errors from bootstrapping.

function [m, res, CI, se, procspeedmean, procspeedCI, param] = fit_gauss_sum(data, param)

if nargin < 2 % defaults
    param.fit_N = 3; % number of Gaussians to fit (2 or 3)
    param.logplot = 0; % fit/plot on logarithmic x scale? (otherwise linear)
    param.bsCI = 1; % use bootstrapping to get confidence intervals?
    param.fixzero = 1; % fix first Gaussian center at 0 nm/s?
end

gauss = @(a,x) a(1)*exp(-(x-a(2)).^2 / (2*a(3))); % Gaussian function

if param.fit_N == 2
    if param.fixzero
        sumgauss = @(a,x) gauss([a(1) 0 a(2)],x) + gauss([a(3) a(4) a(5)],x) + gauss([a(6) a(7) a(8)],x); % Sum of 2 Gaussians, with first fixed at 0
    else
        sumgauss = @(a,x) gauss([a(1) a(2) a(3)],x) + gauss([a(4) a(5) a(6)],x) + gauss([a(7) a(8) a(9)],x); % Sum of 2 Gaussians
    end
elseif param.fit_N == 3
    if param.fixzero
        sumgauss = @(a,x) gauss([a(1) 0 a(2)],x) + gauss([a(3) a(4) a(5)],x) + gauss([a(6) a(7) a(8)],x); % Sum of 3 Gaussians, with first fixed at 0
    else
        sumgauss = @(a,x) gauss([a(1) a(2) a(3)],x) + gauss([a(4) a(5) a(6)],x) + gauss([a(7) a(8) a(9)],x); % Sum of 3 Gaussians
    end
end

if ~param.logplot % linear x axis
    
    histx = 0:2:500;
    [bars, edges] = histcounts(data, histx, 'Normalization', 'Probability');
    
    gaussfun = @(a,x) gauss(a, x);
    fitfun = @(a,x) sumgauss(a, x);
    
    if param.fit_N == 2
        if param.fixzero
            lb = [0 0 0 0 0];
            ub = [1 Inf 1 Inf Inf];
            initval = [0.5 20 0.002 180 8000];
        else
            lb = [0 0 0 0 0 0];
            ub = [1 Inf Inf 1 Inf Inf];
            initval = [0.5 0 20 0.002 180 8000];
        end
    elseif param.fit_N == 3
        if param.fixzero
            lb = [0 0 0 0 0 0 0 0];
            ub = [1 Inf 1 Inf Inf 1 Inf Inf];
            initval = [0.5 20 0.02 27 100 0.002 180 8000];
        else
            lb = [0 0 0 0 0 0 0 0 0];
            ub = [1 Inf Inf 1 Inf Inf 1 Inf Inf];
            initval = [0.5 0 20 0.02 27 100 0.002 180 8000];
        end
    end
    
else % logarithmic x axis (base 10)
    
    histx = log10(0.1):log10(1.5):log10(500);
    [bars, edges] = histcounts(data, 10.^histx, 'Normalization', 'Probability');
    
    gaussfun = @(a,x) gauss(a, log10(x));
    fitfun = @(a,x) sumgauss(a, log10(x));
    
    if param.fit_N == 2
        if param.fixzero
            lb = [0 0 0 0 0];
            ub = [1 Inf 1 Inf Inf];
            initval = [0.1 0.1 0.1 2.3 0.1];
        else
            lb = [0 0 0 0 0 0];
            ub = [1 Inf Inf 1 Inf Inf];
            initval = [0.1 0.1 0.1 0.1 2.3 0.1];
        end
    elseif param.fit_N == 3
        if param.fixzero
            lb = [0 0 0 0 0 0 0 0];
            ub = [10 Inf 10 Inf Inf 10 Inf Inf];
            initval = [0.1 0.1 0.2 1.3 0.1 0.1 2.3 0.1];
        else
            lb = [0 0 0 0 0 0 0 0 0];
            ub = [10 Inf Inf 10 Inf Inf 10 Inf Inf];
            initval = [0.1 0.1 0.1 0.2 1.3 0.1 0.1 2.3 0.1];
        end
    end
    
end

centers = edges(1:end-1) + (edges(2:end)-edges(1:end-1))/2; % want to fit the centers of the bins, not the left edges.

options = optimoptions('lsqcurvefit', 'Display', 'off');
[m, ~, res, ~, ~, ~, J] = lsqcurvefit(fitfun, initval, centers, bars, lb, ub, options);
CI = nlparci(m, res, 'Jacobian', J); % 95% confidence intervals
t = tinv(1-0.05/2, length(bars)-length(m));
se = (CI(:,2)-CI(:,1)) ./ (2*t); % standard errors

tv = 0:0.1:500;
procspeedmean = procspeedmean_fit3gauss(data);
if param.bsCI
    NBOOT = 1000;
    procspeedCI = bootci(NBOOT,@procspeedmean_fit3gauss,data);
else
    NBOOT = [];
    procspeedCI = [];
end

function boot_speed_proc = procspeedmean_fit3gauss(bootdata)
	[bars, edges] = histcounts(bootdata, histx, 'Normalization', 'Probability');
    centers = edges(1:end-1) + (edges(2:end)-edges(1:end-1))/2; % want to fit the centers of the bins, not the left edges.
	[mboot, ~, res, ~, ~, ~, J] = lsqcurvefit(fitfun, initval, centers, bars, lb, ub, options);
	boot_speed_proc = sum(tv.*gaussfun(mboot(3:5),tv))/sum(gaussfun(mboot(3:5),tv));
end

param.hist_xlim = [histx(1) histx(end)];
param.hist_bin_fit = histx(2)-histx(1);
param.fit_fun = fitfun;
param.lower_bounds = lb;
param.upper_bounds = ub;
param.initial_guesses = initval;
param.fit_options = options;
param.N_bootstrap = NBOOT;

end