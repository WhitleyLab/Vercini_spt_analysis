% fit SPT speed data to three independent gaussian functions

logplot = 0; % use logarithmic scale for fitting
figure

data = abs(phmm30_penG.Speed_nm_s_);

[m, res, CI, param] = fit_2gauss(data);

gauss = @(a,x) a(1)*exp(-(x-a(2)).^2 / (2*a(3))); % Gaussian function

gauss1 = @(a,x) gauss([a(1) 0 a(2)],x); % Gaussian function centered at 0
gauss2 = @(a,x) gauss([a(1) 0 a(2)],x) + gauss([a(3) a(4) a(5)],x); % Sum of 2 Gaussians with one centered at 0

if logplot
    [~, edges] = histcounts(log10(data));
    histogram(data, 10.^edges, 'Normalization', 'Probability', 'FaceColor', [0.65 0.65 0.65])
    set(gca, 'xscale', 'log')
    hold on
else
    histx = 0:2:500;
    histogram(data, histx, 'Normalization', 'Probability', 'FaceColor', [0.65 0.65 0.65])
    hold on
end

tv = 0:0.1:500;

plot(tv, gauss1(m(1:2),tv), 'r', 'linew', 2)
plot(tv, gauss(m(3:5),tv), 'c', 'linew', 2)
plot(tv, gauss2(m,tv),'m', 'linew', 2)

xlabel('Speeds (nm/s)')
ylabel('Frequency')
xlim([0 50])
set(gca, 'FontSize', 20)
set(gcf, 'Position', [804 307 950.5 594.5])

param.logplot = logplot;

ud.param = param;
ud.Results = m;
ud.Residuals = res;
ud.CI = CI;
ud.Processive_speed_mean = procspeedmean;
ud.Processive_speed_CI = CI;
ud.Fraction_Processive = proc_frac;
ud.Processive_speed_mean_sanity_check = mean_proc;

set(gcf, 'UserData', ud)
