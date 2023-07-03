% fit SPT speed data to three independent gaussian functions

logfit = 0; % use logarithmic scale for fitting
figure

data = abs(phmm30_penG.Speed_nm_s_);

gam = @(a,x) a(1) * 1/(gamma(a(2))*a(3)^a(2)) * x.^(a(2)-1) .* exp(-x / a(3)); % general gamma function

gamma2 = @(a,x) gam([a(1) a(2) a(3)],x) + gam([a(4) a(5) a(6)],x);
% gamma3 = @(a,x) gam([a(1) a(2) a(3)],x) + gam([a(4) a(5) a(6)],x) + gam([a(7) a(8) a(9)],x);

gamma3 = @(a,x) gam([a(1) 1 a(2)],x) + gam([a(3) a(4) a(5)],x) + gam([a(6) a(7) a(8)],x); % 3 gamma functions where first one is fixed at k=1.

if logfit
    [~, edges] = histcounts(log10(data));
    histogram(data, 10.^edges)
    set(gca, 'xscale', 'log')
    hold on
    
    fitfun = @(a,x) log10(gamma3(a,x));
    lb = [0 0 0 0 0 0 0 0 0];
    ub = [1 Inf Inf 1 Inf Inf 1 Inf Inf];
    initval = [0.5 1 2 0.02 25 100 0.002 180 8000];
else
    histogram(data,0:2:500,'Normalization','Probability')
    hold on
    
    [bars, edges] = histcounts(data,0:2:500,'Normalization','Probability');
    edges = edges(1:end-1) + (edges(2)-edges(1))/2;
    
    fitfun = @(a,x) gamma3(a,x);
%     lb = [0 1 0 0 1 0 0 1 0];
%     ub = [1 Inf Inf 1 Inf Inf 1 Inf Inf];
%     initval = [1 1 10 1 3 10 1 4 40];
    
    lb = [0 0 0 1 0 0 1 0];
    ub = [1 Inf 1 Inf Inf 1 Inf Inf];
    initval = [1 10 1 3 10 1 4 40];
end

% [m, res, J] = nlinfit(edges, bars, gauss3, [0.2 2 0.02 20 4 0.002 200 100]);
% lb = [0 0 0 0 0 0 0 0];
% ub = [1 Inf 1 Inf Inf 1 Inf Inf];
% [m, ~, res, ~, ~, ~, J] = lsqcurvefit(gauss3, [0.5 2 0.02 20 4 0.002 200 100], edges, bars, lb, ub);

% lb = [0 0 0 0 0];
% ub = [1 Inf 1 Inf Inf];
% [m, ~, res, ~, ~, ~, J] = lsqcurvefit(gauss2, [0.5 10 0.04 20 4], edges, bars, lb, ub);
% [~, se] = nlparci2(m, res, 'Jacobian', J);

[m, ~, res, ~, ~, ~, J] = lsqcurvefit(fitfun, initval, edges, bars, lb, ub);
[~, se] = nlparci2(m, res, 'Jacobian', J);

tv = 0:0.1:500;
plot(tv, gam([m(1) 1 m(2)],tv), 'k', 'linew', 1)

% plot(tv, gauss([m(3) 25 100],tv), 'r', 'linew', 1)
% plot(tv, gauss([m(4) 180 8000],tv), 'r', 'linew', 1)

plot(tv, gam(m(3:5),tv), 'r', 'linew', 1)
plot(tv, gam(m(6:8),tv), 'c', 'linew', 1)
plot(tv, fitfun(m,tv),'m', 'linew', 1)