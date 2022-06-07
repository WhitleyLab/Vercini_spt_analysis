function t = add_spt_t_points
%choose set of points (intervals)
%120220 mjc

ud = get(gcf,'UserData');

t = [];

x = 1;
%get points
while ~isempty(x)
    k = waitforbuttonpress;
    if k == 1 %a button was pressed
        x = [];
        return
    else %a mouse press
        x = convx;
        y = convy;
        t = [t x];
        
        % if multiple objects overlapping, figure out which one is closest to the
        % clicks
        h_alive_ind = find(ud.phandles.p_dist~=0); % handles for objects that weren't deleted
        minr=Inf(1,max(h_alive_ind));
        for ii = h_alive_ind
            time_h = get(ud.phandles.p_dist(ii),'xdata');
            dist_h = get(ud.phandles.p_dist(ii),'ydata');
            
            dr = sqrt((time_h-x).^2 + (dist_h-y).^2);

            minr(ii) = min(dr);
        end
        
        % pick closest object
        [~, ind_minr] = min(minr);
        
        %check if list exists (only one could exist)
        if ~isfield(ud.cell.track{ind_minr},'t_pts')
            ud.cell.track{ind_minr}.t_pts = [];
        end
        
        ud.cell.track{ind_minr}.t_pts = [ud.cell.track{ind_minr}.t_pts x];
        set(gcf,'Userdata',ud);
        
        yl_rho = get(ud.axes.ax_rho,'ylim');
        yl_dist = get(ud.axes.ax_dist,'ylim');
        yl_int = get(ud.axes.ax_int,'ylim');
        for ii = 1:length(t)
            plot(ud.axes.ax_rho, t(ii)*[1 1], yl_rho, 'k:', 'linew', 2, 'tag', 'tpts')
            plot(ud.axes.ax_dist, t(ii)*[1 1], yl_dist, 'k:', 'linew', 2, 'tag', 'tpts')
            plot(ud.axes.ax_int, t(ii)*[1 1], yl_int, 'k:', 'linew', 2, 'tag', 'tpts')
        end
        set(gcf,'Userdata',ud);
     end
end

