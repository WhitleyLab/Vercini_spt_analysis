
function resetcenter

ud = get(gcf,'UserData');

k = waitforbuttonpress;
if k == 1 %a button was pressed
    return
else %a mouse press
    [rho, theta] = conv_rho_theta;
end

% polarplot(ud.axes.ax_polar, theta, rho, 'xr')

ud.cell.cell_center_shift_polar = [rho theta];

% original origin in polar coordinates
if isfield(ud.cell,'cell_center_shift_cartesian')
    orix1 = ud.cell.cell_center_shift_cartesian(1);
    oriy1 = ud.cell.cell_center_shift_cartesian(2);
else
    orix1 = 0;
    oriy1 = 0;
end

% new origin in cartesian coordinates
orix2 = rho*cos(theta);
oriy2 = rho*sin(theta);

orix2_2 = orix1 - orix2;
oriy2_2 = oriy1 - oriy2;

ud.cell.cell_center_shift_cartesian = [orix2_2 oriy2_2];

% new origin in polar coordinates
orir2 = sqrt(orix2_2^2 + oriy2_2^2);
orith2 = atan2(oriy2_2, orix2_2);
    
cla(ud.axes.ax_polar)
cla(ud.axes.ax_rho)
cla(ud.axes.ax_dist)

for tt = 1:length(ud.phandles.p_polar)
    
    % original track trajectory in polar coordinates
    r1 = ud.cell.track{tt}.rho_new;
    th1 = ud.cell.track{tt}.theta_new;
    time = ud.cell.track{tt}.time;
    
    % original track trajectory in cartesian coordinates
    x1 = r1.*cos(th1);
    y1 = r1.*sin(th1);

    % new track trajectory in cartesian coordinates
    x2 = x1-orix2;
    y2 = y1-oriy2;
    
    % new track trajectory in polar coordinates
    r2 = sqrt(x2.^2 + y2.^2);
    th2 = atan2(y2, x2);

    ud.cell.track{tt}.rho_new = r2;
    ud.cell.track{tt}.theta_new = th2;
    
    polarplot(ud.axes.ax_polar, th2, r2);
    
    newdist = r2.*th2; % [nm] angular position
    
    % patch up boundary conditions
    dxdt = (newdist(2:end)-newdist(1:end-1)) / ud.param.interval; % first derivative
    ind_plus = find(dxdt>=ud.param.derivative_threshold);
    ind_minus = find(dxdt<=-ud.param.derivative_threshold);
    
    % when track reaches boundary [-pi pi], add or subtract 2*pi*r to keep track continuous
    for jj = 1:length(ind_plus)
        newdist(ind_plus(jj)+1:end) = newdist(ind_plus(jj)+1:end) - 2*pi*r2(ind_plus(jj)+1:end);
    end
    for jj = 1:length(ind_minus)
        newdist(ind_minus(jj)+1:end) = newdist(ind_minus(jj)+1:end) + 2*pi*r2(ind_minus(jj)+1:end);
    end
    
    ud.phandles.p_rho(tt) = plot(ud.axes.ax_rho, time, r2);
    ud.phandles.p_dist(tt) = plot(ud.axes.ax_dist, time, newdist);
    
end

polarplot(ud.axes.ax_polar, orith2, orir2, 'xb')

set(gcf,'userData',ud)

end

function [r, th] = conv_rho_theta
    xy = get(gcf,'currentPoint');
    figpos = get(gcf,'position');
    axpos = get(gca,'position');
    axwidth = axpos(3)*figpos(3); % width of x axis within figure
    axcentx = (axpos(1)+axpos(3)/2)*figpos(3); % zero not at the edge of the axes - in middle because polar axes
    axp = xy(1) - axcentx; % offset x position
    aywidth = axpos(4)*figpos(4); % with of y axis within figure
    aycenty = (axpos(2)+axpos(4)/2)*figpos(4); % zero not at the edge of the axes - in middle because polar axes
    ayp = xy(2) - aycenty; % offset y position
    
%     arp = sqrt(axp^2 + ayp^2); % displacement of click from center
    arp = sqrt((axp/(axwidth/2))^2 + (ayp/(aywidth/2))^2); % displacement of click from center
    meanwidth = mean([axwidth aywidth]);
    
    athp = atan2(-ayp, axp); % [rad] theta. need to reflect about y-axis to match figure
    th = mod(athp, 2*pi); % [rad] convert from [-pi, pi] to [0, 2*pi]
    
    rl = get(gca,'rlim');
%     r = arp/(meanwidth/2) * (rl(2)-rl(1)) + rl(1); % note this needs to be scaled by half the width because center is not at (0,0) of axis position
    r = arp * (rl(2)-rl(1)) + rl(1);
end