function multikey_analysis(src,evnt)
%figure callback function to interactively analyze and interact with a force
%feedback plot
%modified by kw from getpeaks.m (110228 matt comstock) on 220301

ud = get(gcf,'UserData');

if evnt.Character == 'c' % reset center of ring
    resetcenter;
    
elseif evnt.Character == 'a' % run AutoStepFinder on trajectory
    find_steps_spt;
    
elseif evnt.Character == 'b' % run AutoStepFinder on all trajectories (batch version)
    find_steps_spt_batch;
    
elseif evnt.Character == 'd' % run Diffusion Velocity MLE analysis on track
    find_changept_diff_vel;
    
elseif evnt.Character == 't' %add set of t points (for uvrd 2B fret analysis)
    pts = add_spt_t_points;
    
elseif evnt.Character == 'm' %get mean and time for an interval
    get_spt_mean;
    
elseif evnt.Character == 'r' %remove last clicked point if one
    removelastpeak
    
elseif evnt.Character == 's' %save figure
    saveas(gcf,[ud.figsavename '.fig']);
    
elseif evnt.Character == 'v' %calculate velocity between two clicked points
    calc_spt_velocity;

elseif evnt.Character == 'x' % leave this function
    set(gcf,'keypressfcn',[])    
    
end
