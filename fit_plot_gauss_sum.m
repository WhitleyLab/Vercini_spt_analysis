% fit SPT speed data to three independent gaussian functions

param.fit_N = 3;
param.logplot = 0; % use logarithmic scale for fitting
param.bsCI = 0; % bootstrap to get confidence intervals
param.fixzero = 1; % fix the mean of the first Gaussian to be 0 nm/s.

data = abs(phmm30_pc19.Speed_nm_s_);

[fitval, res, CI, se, procspeedmean, procspeedCI, param] = fit_gauss_sum(data, param);

plot_gauss_sum(data, param, fitval, res, CI, se, procspeedmean, procspeedCI)