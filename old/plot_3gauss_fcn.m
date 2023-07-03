% This function plots a histogram of single-molecule speeds and fits a sum
% of three independent Gaussian functions to it. It then calculates the
% expected values of each distribution, with errors from bootstrapping.

% Input logfit = 1 will plot and fit the data on a logarithmic x axis.
% Otherwise it will plot and fit the data on a linear x axis.

function [m, res, CI, procspeedmean, procspeedCI] = plot_3gauss_fcn(data, logfit)

if nargin < 2
    logfit = 0;
end

gauss = @(a,x) a(1)*exp(-(x-a(2)).^2 / (2*a(3))); % general gaussian function

gauss1 = @(a,x) gauss([a(1) 0 a(2)],x); % half-gaussian centered at zero
gauss2 = @(a,x) gauss([a(1) 0 a(2)],x) + gauss([a(3) a(4) a(5)],x);
gauss3 = @(a,x) gauss([a(1) 0 a(2)],x) + gauss([a(3) a(4) a(5)],x) + gauss([a(6) a(7) a(8)],x);

if logfit
    [~, edges] = histcounts(log10(data));
    histogram(data, 10.^edges)
    set(gca, 'xscale', 'log')
    hold on
    
    fitfun = @(a,x) log10(gauss3(a,x));
    lb = [0 0 0 0 0 0 0 0];
    ub = [1 Inf 1 Inf Inf 1 Inf Inf];
    initval = [0.5 2 0.02 25 100 0.002 180 8000];
else
    histx=0:2:500;
    
    histogram(data,histx,'Normalization','Probability')
    hold on

    [bars, edges] = histcounts(data,histx,'Normalization','Probability');
    edges = edges(1:end-1) + (edges(2)-edges(1))/2; % shift bin centers by half a bin width
    
    fitfun = @(a,x) gauss3(a,x);
    lb = [0 0 0 0 0 0 0 0];
    ub = [1 Inf 1 Inf Inf 1 Inf Inf];
    initval = [0.5 20 0.02 27 100 0.002 180 8000];
end

options = optimoptions('lsqcurvefit','Display','off');
[m, ~, res, ~, ~, ~, J] = lsqcurvefit(fitfun, initval, edges, bars, lb, ub, options);
CI = nlparci(m, res, 'Jacobian', J);

tv = 0:0.1:500;
procspeedmean = procspeedmean_fit3gauss(data)
NBOOT=1000;
[procspeedCI,bootstat] = bootci(NBOOT,@procspeedmean_fit3gauss,data);


    function boot_speed_proc=procspeedmean_fit3gauss(bootdata)
        [bars, edges] = histcounts(bootdata,histx,'Normalization','Probability');
        edges = edges(1:end-1) + (edges(2)-edges(1))/2;
        [mboot, ~, res, ~, ~, ~, J] = lsqcurvefit(fitfun, initval, edges, bars, lb, ub,options);
        boot_speed_proc = sum(tv.*gauss(mboot(3:5),tv))/sum(gauss(mboot(3:5),tv));
    end

end