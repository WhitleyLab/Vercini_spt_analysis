% fit SPT speed data to three independent 'Folded Normal' functions

%% Load data

data = abs(cdm30.Speed_nm_s_);

%% Set user parameters

param.logfit = 1; % use logarithmic scale for fitting
param.fit_three = 1; % fit three folded gaussian functions (otherwise only two)
param.fix_gauss1_variance = 0; % fix variance of immobile population, don't fit it

param.FreedmanDiaconis = 0; % use Freedman-Diaconis rule for getting 'optimal' bin widths. otherwise fixed.
bin_width = 2; % [nm/s] widths to use for speed histograms if not using Freedman-Diaconis rule.
% bin_width = 1.26; % [nm/s] widths to use for speed histograms if not using Freedman-Diaconis rule.
bin_width = 1.5; % [nm/s] widths to use for speed histograms if not using Freedman-Diaconis rule.

%% Plot data and fit model to data

fgauss = @(a,x) 1/sqrt(a(2)*2*pi)*exp(-(x-a(1)).^2 / (2*a(2))) + 1/sqrt(a(2)*2*pi)*exp(-(x+a(1)).^2 / (2*a(2))); % folded Gaussian PDF

fgauss1 = @(a,x) fgauss([0 a(1)],x); % folded Gaussian centered at zero
fgauss2 = @(a,x) a(1)*fgauss([0 a(2)],x) + a(3)*fgauss([a(4) a(5)],x);
fgauss3 = @(a,x) a(1)*fgauss([0 a(2)],x) + a(3)*fgauss([a(4) a(5)],x) + a(6)*fgauss([a(7) a(8)],x); % sum of three functions with weights

if param.fix_gauss1_variance % fix variance of immobile population, don't fit it
    fgauss3 = @(a,x) a(1)*fgauss([0 1.5],x) + a(2)*fgauss([a(3) a(4)],x) + a(5)*fgauss([a(6) a(7)],x); % sum of three functions with weights
end

figure

if param.logfit
    if param.FreedmanDiaconis % use Freedman-Diaconis rule to get 'optimal' bin widths
%         [bars, edges] = histcounts(log10(data), 'Normalization', 'Probability', 'BinMethod', 'fd'); % use freedman-diaconis rule for bin sizes
        wid = 2*iqr(log10(data)) / length(log10(data))^(1/3);
    else
        wid = log10(bin_width);
    end
    
    edges = log10(0.1):wid:log10(1000);
    histogram(data, 10.^edges, 'Normalization', 'Probability')
    
%     edges = edges(1:end-1) + (edges(2)-edges(1))/2;
%     histogram(data, 10.^edges, 'Normalization', 'Probability')
    set(gca, 'xscale', 'log')
    hold on

    [bars, edges] = histcounts(log10(data), edges, 'Normalization', 'Probability');
    edges = edges(1:end-1) + (edges(2)-edges(1))/2;
    
    fitfun = @(a,x) fgauss3(a,log10(x));
    edges = 10.^edges;
    lb = [0 0 0 0 0 0 0 0];
    ub = [10 Inf 10 Inf Inf 10 Inf Inf];
    initval = [0.1 0.1 0.2 1.3 0.1 0.1 2.3 0.1];
    
    fgauss1 = @(a,x) fgauss1(a,log10(x));    
    fgauss = @(a,x) fgauss(a,log10(x));
    
    xlims = [0.1 1e3];
else
    if param.FreedmanDiaconis% use Freedman-Diaconis rule to get 'optimal' bin widths
        wid = 2*iqr(data) / length(data)^(1/3);
    else
        wid = bin_width;
    end
    
    histogram(data,0:wid:500,'Normalization','Probability')
    hold on
    
    [bars, edges] = histcounts(data,0:wid:500,'Normalization','Probability');
    edges = edges(1:end-1) + (edges(2)-edges(1))/2;
    
    if param.fit_three
        fitfun = @(a,x) fgauss3(a,x);
        if param.fix_gauss1_variance
            lb = [0 0 0 0 0 0 0];
            ub = [10 10 Inf Inf 10 Inf Inf];
            initval = [1 1 15 50 1 150 8000];
        else
            lb = [0 0 0 0 0 0 0 0];
            ub = [10 Inf 10 Inf Inf 10 Inf Inf];
            initval = [1 1.5 1 15 50 1 150 8000];
        end
    else
        fitfun = @(a,x) fgauss2(a,x);
        if param.fix_gauss1_variance
            lb = [0 0 0 0 0];
            ub = [100 Inf 10 Inf Inf];
            initval = [1 1.5 1 150 8000];
        else
            lb = [0 0 0 0];
            ub = [100 10 Inf Inf];
            initval = [1 1 150 8000];
        end
    end
    xlims = [0 100];

end

[m, ~, res, ~, ~, ~, J] = lsqcurvefit(fitfun, initval, edges, bars, lb, ub);
[~, se] = nlparci2(m, res, 'Jacobian', J);

Results = [];
if param.fix_gauss1_variance
    Results.Amplitude1 = m(1);
    Results.Amplitude1_err = se(1);
    Results.Amplitude2 = m(2);
    Results.Amplitude2_err = se(2);
    Results.Center2 = m(3);
    Results.Center2_err = se(3);
    Results.Variance2 = m(4);
    Results.Variance2_err = se(4);
    if param.fit_three
        Results.Amplitude3 = m(5);
        Results.Amplitude3_err = se(5);
        Results.Center3 = m(6);
        Results.Center3_err = se(6);
        Results.Variance3 = m(7);
        Results.Variance3_err = se(7);
    end
else
    Results.Amplitude1 = m(1);
    Results.Amplitude1_err = se(1);
    Results.Variance1 = m(2);
    Results.Variance1_err = se(2);
    Results.Amplitude2 = m(3);
    Results.Amplitude2_err = se(3);
    Results.Center2 = m(4);
    Results.Center2_err = se(4);
    Results.Variance2 = m(5);
    Results.Variance2_err = se(5);
    if param.fit_three
        Results.Amplitude3 = m(6);
        Results.Amplitude3_err = se(6);
        Results.Center3 = m(7);
        Results.Center3_err = se(7);
        Results.Variance3 = m(8);
        Results.Variance3_err = se(8);
    end
end

%% Plot model

if param.logfit
    plotx_bin = wid/10;
    plotx = edges(1):plotx_bin:edges(end);
else
    plotx = 0:0.1:500;
end

if param.fix_gauss1_variance
    plot(plotx, m(1)*fgauss1(1.5,plotx), 'k', 'linew', 1)
    
    plot(plotx, m(2)*fgauss(m(3:4),plotx), 'r', 'linew', 1)
    if param.fit_three
        plot(plotx, m(5)*fgauss(m(6:7),plotx), 'c', 'linew', 1)
    end
else
    plot(plotx, m(1)*fgauss1(m(2),plotx), 'k', 'linew', 1)
    
    plot(plotx, m(3)*fgauss(m(4:5),plotx), 'r', 'linew', 1)
    if param.fit_three
        plot(plotx, m(6)*fgauss(m(7:8),plotx), 'c', 'linew', 1)
    end
end
plot(plotx, fitfun(m,plotx),'m', 'linew', 1)

% get fractions of each population. Need to divide by bin width, because
% histograms were normalized by a sum, whereas we want the area under each
% curve (which is an integral, not a sum).
Results.Fraction1 = Results.Amplitude1 / wid;
Results.Fraction2 = Results.Amplitude2 / wid;
if param.fit_three
    Results.Fraction3 = Results.Amplitude3 / wid;
end

Results.Fraction1_err = se(1) / wid;
Results.Fraction2_err = se(3) / wid;
if param.fit_three
    Results.Fraction3_err = se(6) / wid;
end

xlim(xlims)
xlabel('Speed (nm/s)')
ylabel('Frequency')
title('Speed distribution for HaloTag-PBP2B')

%% Save parameters and data into figure.

param.bin_widths = wid; % only care what the final one is.
param.fit_function = fitfun;
param.lower_bounds = lb;
param.upper_bounds = ub;
param.initial_fit_params = initval;

ud.param = param;
ud.data = data;
ud.results = Results;

set(gcf, 'UserData', ud)