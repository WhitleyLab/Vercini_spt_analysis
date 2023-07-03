function find_changept_diff_vel
%click near a trajectory and run diffusion-velocity changepoint detection on it.

ud = get(gcf,'UserData');

k = waitforbuttonpress;
if k == 1 %a button was pressed
    return
else %a mouse press
    x = convx;
    y = convy;
end

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

time = get(ud.phandles.p_dist(ind_minr),'xdata');
dist = get(ud.phandles.p_dist(ind_minr),'ydata');

dy = dist(2:end)-dist(1:end-1);

codepath = 'C:\Users\nkw81\Documents\GitHub\DiffusionVelocityChangePoint\parallel-DV\';

% write dy coordinates to temporary txt file for diffusion-velocity changepoint code to use
fileid = fopen([ud.figsavename '.txt'], 'w');
fprintf(fileid, '%f\r\n', dy./1000); % convert to um first, just in case that's important
fclose(fileid);

alpha = 0.05;
system(['mpiexec -np 2 ' codepath 'Parallel-DV "' ud.figsavename '.txt" ' num2str(ud.param.interval) ' ' num2str(alpha)]);

% get data from new file
outfile = fopen([ud.figsavename '.txt.cp'], 'r');
chgpts = fscanf(outfile, '%d');
fclose(outfile);

ud.cell.track{ind_minr}.diff_vel_changepoints = chgpts;

% delete temporary files
delete([ud.figsavename '.txt'])
delete([ud.figsavename '.txt.cp'])

set(gcf,'userData',ud)
