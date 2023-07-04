function plot_spt_tpts
%plot the time pt intervals for fret
ud = get(gcf,'UserData');

%clear previous plots
clear_plots;

t = ud.t_pts;

%now plot
%plot on FRET plot, markers first
if isfield(ud.param,'fret_axh')
    axes(ud.param.fret_axh);
    hold all
    
    % if ~isempty(t)
    %     ud.param.t_pts_plot_h = plot(t,t*0+0.5,'k+','linew',2,'tag','fpts')
    % end
    
    %dashed lines
    %on fret plot
    for i = 1:length(t)
        plot(t(i)*[1 1],[-1 2],'k:','linew',2,'tag','fpts')
    end
end

%on trap plot
axes(ud.phandles.axesh)
yl = get(gca,'ylim');
hold all
for i = 1:length(t)
    plot(t(i)*[1 1],yl,'k:','linew',2,'tag','fpts')
end

%on fluorescence plot
if isfield(ud.param,'fl_axh')
    axes(ud.param.fl_axh)
    yl = get(gca,'ylim');
    hold all
    for i = 1:length(t)
        plot(t(i)*[1 1],yl,'k:','linew',2,'tag','fpts')
    end
end

set(gcf,'Userdata',ud);