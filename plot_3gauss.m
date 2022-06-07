% fit SPT speed data to three independent gaussian functions

data = abs(phmm37.Speed_nm_s_);

% gauss = @(a,x) heaviside(x) .* (a(1)*exp(-(x-a(2)).^2 / (2*a(3)))); % general gaussian function, bounded at zero
gauss = @(a,x) a(1)*exp(-(x-a(2)).^2 / (2*a(3))); % general gaussian function, bounded at zero

gauss1 = @(a,x) gauss([a(1) 0 a(2)],x); % half-gaussian centered at zero
gauss2 = @(a,x) gauss([a(1) 0 a(2)],x) + gauss([a(3) a(4) a(5)],x);
gauss3 = @(a,x) gauss([a(1) 0 a(2)],x) + gauss([a(3) a(4) a(5)],x) + gauss([a(6) a(7) a(8)],x);

[bars, edges] = histcounts(data,0:1:500,'Normalization','Probability');
edges = edges(1:end-1) + (edges(2)-edges(1))/2;

% [m, res, J] = nlinfit(edges, bars, gauss3, [0.2 2 0.02 20 4 0.002 200 100]);
% lb = [0 0 0 0 0 0 0 0];
% ub = [1 Inf 1 Inf Inf 1 Inf Inf];
% [m, ~, res, ~, ~, ~, J] = lsqcurvefit(gauss3, [0.5 2 0.02 20 4 0.002 200 100], edges, bars, lb, ub);

lb = [0 0 0 0 0];
ub = [1 Inf 1 Inf Inf];
[m, ~, res, ~, ~, ~, J] = lsqcurvefit(gauss2, [0.5 10 0.04 20 4], edges, bars, lb, ub);
[~, se] = nlparci2(m, res, 'Jacobian', J);

tv = 0:0.1:500;
plot(tv, gauss1(m(1:2),tv), 'k', 'linew', 1)
plot(tv, gauss(m(3:5),tv), 'r', 'linew', 1)
% plot(tv, gauss(m(6:8),tv), 'c', 'linew', 1)