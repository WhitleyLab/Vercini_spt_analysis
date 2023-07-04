% fit SPT speed data to three independent gaussian functions

logfit = 0; % use logarithmic scale for fitting
figure

data = abs(t111a.Speed_nm_s_);

gauss = @(a,x) heaviside(x) .* (a(1)*exp(-(x-a(2)).^2 / (2*a(3)))); % truncated gaussian function (bounded at zero)
% gauss = @(a,x) a(1)*exp(-(x-a(2)).^2 / (2*a(3))); % general gaussian function

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
    histogram(data,0:2:500,'Normalization','Probability')
    hold on
    
    [bars, edges] = histcounts(data,0:2:500,'Normalization','Probability');
    edges = edges(1:end-1) + (edges(2)-edges(1))/2;
    
    fitfun = @(a,x) gauss3(a,x);
    lb = [0 0 0 0 0 0 0 0];
    ub = [1 Inf 1 Inf Inf 1 Inf Inf];
    initval = [0.5 20 0.02 27 100 0.002 180 8000];
    
%     fitfun = @(a,x) gauss3([a(1) a(2) a(3) 25 100 a(4) 180 8000],x);
%     lb = [0 0 0 0];
%     ub = [1 Inf 1 1];
%     initval = [0.5 2 0.02 0.002];
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
plot(tv, gauss1(m(1:2),tv), 'k', 'linew', 1)

% plot(tv, gauss([m(3) 25 100],tv), 'r', 'linew', 1)
% plot(tv, gauss([m(4) 180 8000],tv), 'r', 'linew', 1)

plot(tv, gauss(m(3:5),tv), 'r', 'linew', 1)
plot(tv, gauss(m(6:8),tv), 'c', 'linew', 1)
plot(tv, fitfun(m,tv),'m', 'linew', 1)
xlabel('Speeds (nm/s)')
ylabel('Frequency')
xlim([0 50])