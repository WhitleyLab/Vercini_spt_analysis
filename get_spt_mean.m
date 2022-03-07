function get_spt_mean
%click two points and calculate the mean value between them and the time
%interval
%modified by kw from getmean.m (120626 matt comstock) on 220301

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
minr1=[]; minr2=[];
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
rho = get(ud.phandles.p_rho(ind_minr),'ydata');

[~, it1] = min(abs(time - x1));
[~, it2] = min(abs(time - x2));

t1 = time(it1);
t2 = time(it2);

dwell = t2 - t1;

meandist = mean(dist(it1:it2));
meanrho = mean(rho(it1:it2));

disp(['<position> = ' num2str(meandist,'%0.1f') ' bp, ' num2str(dwell,'%0.1f') ' s'])
disp(['<radius> = ' num2str(meanrho,'%0.1f') ' nm, ' num2str(dwell,'%0.1f') ' s'])

plot(ud.axes.ax_dist, [t1 t2], meandist*[1 1], 'b:', 'tag', 'meandist', 'linew', 2)
plot(ud.axes.ax_rho, [t1 t2], meanrho*[1 1], 'b:', 'tag', 'meanrho', 'linew', 2)

text(ud.axes.ax_dist, t1+(t2-t1)/2,meandist+5,{[num2str(meandist,'%0.1f') ' nm'] [num2str(dwell,'%0.1f') ' s']},'tag','meandist','color','b')
text(ud.axes.ax_rho, t1+(t2-t1)/2,meanrho+5,{[num2str(meanrho,'%0.1f') ' nm'] [num2str(dwell,'%0.1f') ' s']},'tag','meanrho','color','b')

end
