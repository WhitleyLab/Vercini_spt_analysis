% Author: Kevin Whitley
% Date created: 230703

% This function plots a histogram of speeds (in [nm/s]) along with the
% associated fits to a sum of Gaussians. Input parameters and results of
% fitting are saved into the figure as UserData variable ud.
% 
% INPUTS:
% - data: an array of speeds in [nm/s]
% - param: a structure variable with (at least) three fields:
%       - fit_N: The number of independent Gaussian functions you want to
%       fit (either 2 or 3).
%       - logplot: Whether you want to fit and plot the data on a
%       logarithmic scale (speeds as log base 10). Set to 1 for logarithmic
%       plots, or 0 for linear.
%       - fixzero: Whether you want to fix the center of the first Gaussian
%       to be 0 nm/s. Otherwise this will be fitted. Set to 1 to fix this
%       value, or 0 to let it float.
% - fitval: The fitted values from a sum of Gaussian fits. Depending on
% what options you chose in parameters (fit_N and fixzero), this will be
% between 5 and 9 values.
% - res: Residuals from the fits.
% - CI: 95% confidence intervals of fitted parameters.
% - se: Standard error of the mean of fitted parameters.
% - procspeedmean: Mean of the processive (middle) population, calculated
% as the first moment of the fitted Gaussian.
% - procspeedCI: 95% confidence intervals of the processive (middle)
% population, calculated by bootstrapping.

function plot_gauss_sum(data, param, fitval, res, CI, se, procspeedmean, procspeedCI)

gauss = @(a,x) a(1)*exp(-(x-a(2)).^2 / (2*a(3))); % Gaussian function

if param.fit_N ==2
    if param.fixzero
        gauss1 = @(a,x) gauss([a(1) 0 a(2)],x); % Gaussian function centered at 0
        sumgauss = @(a,x) gauss([a(1) 0 a(2)],x) + gauss([a(3) a(4) a(5)],x); % Sum of 2 Gaussians with one centered at 0
    else
        gauss1 = @(a,x) gauss([a(1) a(2) a(3)],x);
        sumgauss = @(a,x) gauss([a(1) a(2) a(3)],x) + gauss([a(4) a(5) a(6)],x); % Sum of 2 Gaussians
    end
elseif param.fit_N ==3
    if param.fixzero
        gauss1 = @(a,x) gauss([a(1) 0 a(2)],x); % Gaussian function centered at 0
        sumgauss = @(a,x) gauss([a(1) 0 a(2)],x) + gauss([a(3) a(4) a(5)],x) + gauss([a(6) a(7) a(8)],x); % Sum of 3 Gaussians with one centered at 0
    else
        gauss1 = @(a,x) gauss([a(1) a(2) a(3)],x);
        sumgauss = @(a,x) gauss([a(1) a(2) a(3)],x) + gauss([a(4) a(5) a(6)],x) + gauss([a(7) a(8) a(9)],x); % Sum of 3 Gaussians
    end
end

if param.logplot
    histx = log10(0.1):log10(1.5):log10(1000);
    
    histogram(data, 10.^histx, 'Normalization', 'Probability', 'FaceColor', [0.65 0.65 0.65])
    set(gca, 'xscale', 'log')
    hold on
    
    tv = 0.1:0.1:1000;
    gauss1 = @(a,x) gauss1(a, log10(x));
    sumgauss = @(a,x) sumgauss(a, log10(x));
    gauss = @(a,x) gauss(a, log10(x));
    
    if param.fit_N == 2
        if param.fixzero
            plot(tv, gauss1(fitval(1:2), tv), 'r', 'linew', 2)
            plot(tv, gauss(fitval(3:5), tv), 'c', 'linew', 2)
            plot(tv, sumgauss(fitval, tv),'m', 'linew', 2)
        else
            plot(tv, gauss1(fitval(1:3), tv), 'r', 'linew', 2)
            plot(tv, gauss(fitval(4:6), tv), 'c', 'linew', 2)
            plot(tv, sumgauss(fitval, tv),'m', 'linew', 2)
        end
    elseif param.fit_N ==3
        if param.fixzero
            plot(tv, gauss1(fitval(1:2), tv), 'r', 'linew', 2)
            plot(tv, gauss(fitval(3:5), tv), 'Color', [0 0.45 0.74], 'linew', 2)
            plot(tv, gauss(fitval(6:8), tv), 'c', 'linew', 2)
            plot(tv, sumgauss(fitval, tv),'m', 'linew', 2)
        else
            plot(tv, gauss1(fitval(1:3), tv), 'r', 'linew', 2)
            plot(tv, gauss(fitval(4:6), tv), 'Color', [0 0.45 0.74], 'linew', 2)
            plot(tv, gauss(fitval(7:9), tv), 'c', 'linew', 2)
            plot(tv, sumgauss(fitval, tv),'m', 'linew', 2)
        end
    end
    
    xlim([1e-1 1e3])
else
    histx = 0:2:500;
    histogram(data, histx, 'Normalization', 'Probability', 'FaceColor', [0.65 0.65 0.65])
    hold on
    
    tv = 0:0.1:500;

    if param.fit_N == 2
        if param.fixzero
            plot(tv, gauss1(fitval(1:2), tv), 'r', 'linew', 2)
            plot(tv, gauss(fitval(3:5), tv), 'Color', [0 0.45 0.74], 'linew', 2)
            plot(tv, sumgauss(fitval, tv),'m', 'linew', 2)
        else
            plot(tv, gauss1(fitval(1:3), tv), 'r', 'linew', 2)
            plot(tv, gauss(fitval(4:6), tv), 'Color', [0 0.45 0.74], 'linew', 2)
            plot(tv, sumgauss(fitval, tv),'m', 'linew', 2)
        end
    elseif param.fit_N ==3
        if param.fixzero
            plot(tv, gauss1(fitval(1:2), tv), 'r', 'linew', 2)
            plot(tv, gauss(fitval(3:5), tv), 'Color', [0 0.45 0.74], 'linew', 2)
            plot(tv, gauss(fitval(6:8), tv), 'c', 'linew', 2)
            plot(tv, sumgauss(fitval, tv),'m', 'linew', 2)
        else
            plot(tv, gauss1(fitval(1:3), tv), 'r', 'linew', 2)
            plot(tv, gauss(fitval(4:6), tv), 'Color', [0 0.45 0.74], 'linew', 2)
            plot(tv, gauss(fitval(7:9), tv), 'c', 'linew', 2)
            plot(tv, sumgauss(fitval, tv),'m', 'linew', 2)
        end
    end
    
    xlim([0 50])
end

xlabel('Speeds (nm/s)')
ylabel('Frequency')

set(gca, 'FontSize', 20)
set(gcf, 'Position', [804 307 950.5 594.5])

if param.fixzero
    area1 = trapz(tv, gauss1(fitval(1:2),tv));%static
    area2 = trapz(tv, gauss(fitval(3:5),tv));%processive
    mean_proc = sum(tv.*gauss(fitval(3:5),tv))/sum(gauss(fitval(3:5),tv));
else
    area1 = trapz(tv, gauss1(fitval(1:3),tv));%static
    area2 = trapz(tv, gauss(fitval(4:6),tv));%processive
    mean_proc = sum(tv.*gauss(fitval(4:6),tv))/sum(gauss(fitval(4:6),tv));
end

proc_frac = area2/(area1+area2);%ignoring the diffusive population

ud = [];
ud.param = param;
ud.Results = fitval;
ud.Residuals = res;
ud.CI = CI;
ud.SEM = se;
ud.Processive_speed_mean = procspeedmean;
ud.Processive_speed_CI = procspeedCI;
ud.Fraction_Processive = proc_frac;
ud.Processive_speed_mean_sanity_check = mean_proc;

set(gcf, 'UserData', ud)
