% plot line segment from autostepfinder parameters

function plot_line_segment_ASF(t, dist, interval, results, ax_hand, wdw)

lineEqn = @(a,x) a(1) + a(2)*x;

speeds = [results.LevelBefore; results.LevelAfter(end)] / interval; % [nm/s]
wind_shift = floor(wdw/2); % derivative used sliding window, so 'true' transition times are slightly shifted behind
time_shift = t(1) + wind_shift - 0.5; % [s] subtracting 0.5 because the transitions occurs between frames
times = [t(1); results.TimeStep+time_shift; t(end)];

v_seg=[]; b_seg=[]; Pos=[]; Time=[];
for ii = 1:length(times)-1
    
    t_seg = times(ii):interval/2:times(ii+1);
    v_seg(ii) = speeds(ii);
    
    % get y-intercept for each line segment
    if ii>1
        b_seg(ii) = (v_seg(ii-1)-v_seg(ii))*Time(end) + b_seg(ii-1);
    else % first one is tricky - need to do a quick linear regression on this segment to find an appropriate y-intercept
        t_seg_all = t(t>=times(ii) & t<=times(ii+1));
        X = [t_seg_all' ones(length(t_seg_all),1)];
        d_seg = dist(1:length(t_seg_all));
        b = nlinfit(X(:,1), d_seg', @(a,x)lineEqn([a v_seg(ii)],x), dist(1));
        b_seg(ii) = b;
    end
    
    Time = [Time t_seg];
    Pos = [Pos lineEqn([b_seg(ii) v_seg(ii)], t_seg)];

end

plot(ax_hand, Time, Pos, 'k')