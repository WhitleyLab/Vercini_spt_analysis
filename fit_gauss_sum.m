% Authors: Kevin Whitley and Séamus Holden
% Date created: 230628
% 
% This function fits a sum of independent Gaussian functions to an array of
% speed data (in [nm/s]). It then calculates the expected values of the
% processive (middle) population. Confidence intervals are calculatedfrom
% bootstrapping.
% 
% INPUTS:
% - data: an array of speeds in [nm/s]
% - param: a structure variable with four fields:
%       - fit_N: The number of independent Gaussian functions you want to
%       fit (either 2 or 3).
%       - logplot: Whether you want to fit and plot the data on a
%       logarithmic scale (speeds as log base 10). Set to 1 for logarithmic
%       plots, or 0 for linear.
%       - bsCI: Whether you want to get 95% confidence intervals for the
%       processive population via bootstrapping. Set to 1 for yes and 0 for
%       no. (Note that bootstrapping takes a fairly long time).
%       - fixzero: Whether you want to fix the center of the first Gaussian
%       to be 0 nm/s. Otherwise this will be fitted. Set to 1 to fix this
%       value, or 0 to let it float.
% 
% OUTPUTS:
% - fitval: The fitted values from the sum of Gaussian fits. Depending on
% what options you chose in parameters (fit_N and fixzero), this will be
% between 5 and 9 values.
% - res: Residuals from the fits.
% - CI: 95% confidence intervals of fitted parameters.
% - se: Standard error of the mean of fitted parameters.
% - procspeedmean: Mean of the processive (middle) population, calculated
% as the first moment of the fitted Gaussian.
% - procspeedCI: 95% confidence intervals of the processive (middle)
% population, calculated by bootstrapping. If bsCI was set to 0, this will
% return an empty array.
% - param: A structure variable containing the inputted parameters as well
% as any that were used here. This can be saved later with any plots so you
% remember what parameters produced it.

function [fitval, res, CI, se, procspeedmean, procspeedCI, param] = fit_gauss_sum(data, param)

if nargin < 2 % defaults
    param.fit_N = 3; % number of Gaussians to fit (2 or 3)
    param.logplot = 0; % fit/plot on logarithmic x scale? (otherwise linear)
    param.bsCI = 1; % use bootstrapping to get confidence intervals?
    param.fixzero = 1; % fix first Gaussian center at 0 nm/s?
end

%% Set fit function(s) depending on inputted parameters

gauss = @(a,x) a(1)*exp(-(x-a(2)).^2 / (2*a(3))); % Gaussian function

if param.fit_N == 2
    if param.fixzero
        sumgauss = @(a,x) gauss([a(1) 0 a(2)],x) + gauss([a(3) a(4) a(5)],x); % Sum of 2 Gaussians, with first fixed at 0
    else
        sumgauss = @(a,x) gauss([a(1) a(2) a(3)],x) + gauss([a(4) a(5) a(6)],x); % Sum of 2 Gaussians
    end
elseif param.fit_N == 3
    if param.fixzero
        sumgauss = @(a,x) gauss([a(1) 0 a(2)],x) + gauss([a(3) a(4) a(5)],x) + gauss([a(6) a(7) a(8)],x); % Sum of 3 Gaussians, with first fixed at 0
    else
        sumgauss = @(a,x) gauss([a(1) a(2) a(3)],x) + gauss([a(4) a(5) a(6)],x) + gauss([a(7) a(8) a(9)],x); % Sum of 3 Gaussians
    end
end

%% Prepare data for fitting, set fit parameters depending on inputted parameters

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
            ub = [10 Inf 10 Inf Inf];
            initval = [0.1 0.1 0.01 2.2 0.1];
        else
            lb = [0 0 0 0 0 0];
            ub = [10 Inf Inf 10 Inf Inf];
            initval = [0.1 0.1 0.1 0.01 2.2 0.1];
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

%% Fit data

options = optimoptions('lsqcurvefit', 'Display', 'off');
[fitval, ~, res, ~, ~, ~, J] = lsqcurvefit(fitfun, initval, centers, bars, lb, ub, options);
CI = nlparci(fitval, res, 'Jacobian', J); % 95% confidence intervals
t = tinv(1-0.05/2, length(bars)-length(fitval));
se = (CI(:,2)-CI(:,1)) ./ (2*t); % standard errors

%% Bootstrap to get 95% confidence intervals
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
    
    if param.fit_N ==3 % only really need to do this if fitting 3 gaussians. no processive population otherwise.
        if param.fixzero
            boot_speed_proc = sum(tv.*gaussfun(mboot(3:5),tv))/sum(gaussfun(mboot(3:5),tv));
        else
            boot_speed_proc = sum(tv.*gaussfun(mboot(4:6),tv))/sum(gaussfun(mboot(4:6),tv));
        end
    else
        boot_speed_proc = [];
    end
end

%% Add all relevant parameters to param for output

param.hist_xlim = [histx(1) histx(end)];
param.hist_bin_fit = histx(2)-histx(1);
param.fit_fun = fitfun;
param.lower_bounds = lb;
param.upper_bounds = ub;
param.initial_guesses = initval;
param.fit_options = options;
param.N_bootstrap = NBOOT;

end