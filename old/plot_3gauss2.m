% fit SPT speed data to three independent gaussian functions

param.fit_N = 3;
param.logplot = 0; % use logarithmic scale for fitting
param.bsCI = 0; % bootstrap to get confidence intervals
param.fixzero = 0; % fix the mean of the first Gaussian to be 0 nm/s.

figure

data = abs(phmm37.Speed_nm_s_);

[m, res, CI, se, procspeedmean, procspeedCI, param] = fit_3gauss(data, param);

gauss = @(a,x) a(1)*exp(-(x-a(2)).^2 / (2*a(3))); % Gaussian function

if param.fixzero
    gauss1 = @(a,x) gauss([a(1) 0 a(2)],x); % Gaussian function centered at 0
    gauss2 = @(a,x) gauss([a(1) 0 a(2)],x) + gauss([a(3) a(4) a(5)],x);
    gauss3 = @(a,x) gauss([a(1) 0 a(2)],x) + gauss([a(3) a(4) a(5)],x) + gauss([a(6) a(7) a(8)],x); % Sum of 3 Gaussians with one centered at 0
else
    gauss1 = @(a,x) gauss([a(1) a(2) a(3)],x);
    gauss2 = @(a,x) gauss([a(1) a(2) a(3)],x) + gauss([a(4) a(5) a(6)],x);
    gauss3 = @(a,x) gauss([a(1) a(2) a(3)],x) + gauss([a(4) a(5) a(6)],x) + gauss([a(7) a(8) a(9)],x);
end

if param.logplot
    histx = log10(0.1):log10(1.5):log10(1000);

    histogram(data, 10.^histx, 'Normalization', 'Probability', 'FaceColor', [0.65 0.65 0.65])
    set(gca, 'xscale', 'log')
    hold on

    tv = 0.1:0.1:1000;
    gauss1 = @(a,x) gauss1(a, log10(x));
    gauss3 = @(a,x) gauss3(a, log10(x));
    gauss = @(a,x) gauss(a, log10(x));

    if param.fixzero
        plot(tv, gauss1(m(1:2), tv), 'r', 'linew', 2)
        plot(tv, gauss(m(3:5), tv), 'Color', [0 0.45 0.74], 'linew', 2)
        plot(tv, gauss(m(6:8), tv), 'c', 'linew', 2)
        plot(tv, gauss3(m, tv),'m', 'linew', 2)
    else
        plot(tv, gauss1(m(1:3), tv), 'r', 'linew', 2)
        plot(tv, gauss(m(4:6), tv), 'Color', [0 0.45 0.74], 'linew', 2)
        plot(tv, gauss(m(7:9), tv), 'c', 'linew', 2)
        plot(tv, gauss3(m, tv),'m', 'linew', 2)
    end
    
    xlim([1e-1 1e3])
else
    histx = 0:2:500;
    histogram(data, histx, 'Normalization', 'Probability', 'FaceColor', [0.65 0.65 0.65])
    hold on
    
    tv = 0:0.1:500;

    if param.fixzero
        plot(tv, gauss1(m(1:2), tv), 'r', 'linew', 2)
        plot(tv, gauss(m(3:5), tv), 'Color', [0 0.45 0.74], 'linew', 2)
        plot(tv, gauss(m(6:8), tv), 'c', 'linew', 2)
        plot(tv, gauss3(m, tv),'m', 'linew', 2)
    else
        plot(tv, gauss1(m(1:3), tv), 'r', 'linew', 2)
        plot(tv, gauss(m(4:6), tv), 'Color', [0 0.45 0.74], 'linew', 2)
        plot(tv, gauss(m(7:9), tv), 'c', 'linew', 2)
        plot(tv, gauss3(m, tv),'m', 'linew', 2)
    end
    
    xlim([0 50])
end

xlabel('Speeds (nm/s)')
ylabel('Frequency')

set(gca, 'FontSize', 20)
set(gcf, 'Position', [804 307 950.5 594.5])

if param.fixzero
    area1 = trapz(tv, gauss1(m(1:2),tv));%static
    area2 = trapz(tv, gauss(m(3:5),tv));%processive
    mean_proc = sum(tv.*gauss(m(3:5),tv))/sum(gauss(m(3:5),tv));
else
    area1 = trapz(tv, gauss1(m(1:3),tv));%static
    area2 = trapz(tv, gauss(m(4:6),tv));%processive
    mean_proc = sum(tv.*gauss(m(4:6),tv))/sum(gauss(m(4:6),tv));
end

proc_frac = area2/(area1+area2);%ignoring the diffusive population

ud = [];
ud.param = param;
ud.Results = m;
ud.Residuals = res;
ud.CI = CI;
ud.SEM = se;
ud.Processive_speed_mean = procspeedmean;
ud.Processive_speed_CI = procspeedCI;
ud.Fraction_Processive = proc_frac;
ud.Processive_speed_mean_sanity_check = mean_proc;

set(gcf, 'UserData', ud)
