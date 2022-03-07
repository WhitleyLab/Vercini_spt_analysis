function calc_spt_velocity
%click two points and calculate the unwinding velocity between the points
%120626 mjc

ud = get(gcf,'UserData');

%get two t points
k = waitforbuttonpress;
if k == 1 %a button was pressed
    return
else %a mouse press
    x1 = convx;
    y1 = convy;
end

k = waitforbuttonpress;
if k == 1 %a button was pressed
    return
else %a mouse press
    x2 = convx;
    y2 = convy;
end

% if multiple objects overlapping, figure out which one is closest to the
% clicks
h_alive_ind = find(ud.phandles.p_dist~=0); % handles for objects that weren't deleted
minr1=Inf(1,max(h_alive_ind)); minr2=Inf(1,max(h_alive_ind));
for ii = h_alive_ind
    time_h = get(ud.phandles.p_dist(ii),'xdata');
    dist_h = get(ud.phandles.p_dist(ii),'ydata');

    dr1 = sqrt((time_h-x1).^2 + (dist_h-y1).^2);
    dr2 = sqrt((time_h-x2).^2 + (dist_h-y2).^2);
    
    minr1(ii) = min(dr1);
    minr2(ii) = min(dr2);
end

% pick closest object
[~, ind_minr] = min(minr1 + minr2);

time = get(ud.phandles.p_dist(ind_minr),'xdata');
dist = get(ud.phandles.p_dist(ind_minr),'ydata');

[~, it1] = min(abs(time - x1));
[~, it2] = min(abs(time - x2));

t1 = time(it1);
t2 = time(it2);

d1 = dist(it1);
d2 = dist(it2);

%simple velocity calculation using endpoints
v_endpoints = (d2-d1)/(t2-t1);
display(['v_endpoints = ' num2str(v_endpoints,'%0.4f')])
plot([t1 t2],[d1 d2],'b:','tag','vel','linew',2)

%fit a line to the interval to get the velocity
x = time(it1:it2);
y = dist(it1:it2);

f = fit(x',y','poly1');
v_fit = f.p1;
display(['v_fit = ' num2str(v_fit,'%0.2f')])
plot(ud.axes.ax_dist, [t1 t2], f.p1*[t1 t2]+f.p2, 'k:', 'tag', 'vel', 'linew',2)

%keep text at top
if v_fit < 0
    txty = d1 + 3;
else
    txty = d2 + 3;
end
text(ud.axes.ax_dist, t1, txty, {[num2str(v_endpoints,'%0.1f') ' nm/s'] [num2str(v_fit,'%0.1f') ' nm/s (fit)']}, 'tag', 'vel', 'color', 'b')

%add fit data to list saved to figure

%check if list exists (only one could exist)
if ~isfield(ud.cell.track{ind_minr},'speeds')
    ud.cell.track{ind_minr}.speeds = [];
end

allv_res = [d1 d2 t1 t2 v_fit];
ud.cell.track{ind_minr}.speeds = [ud.cell.track{ind_minr}.speeds; allv_res];

set(gcf,'userData',ud)

end