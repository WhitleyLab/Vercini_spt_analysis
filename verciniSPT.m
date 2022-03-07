
function verciniSPT(directory, param)

if nargin == 0
    param.Date = '211215';
    param.file = '6';
    directory = ['\\campus\rdw\FMS CBCB\nsh167\Shared\data\Whitley-Kevin\211215_sh147_phmm_37c_vercini\' param.Date '_' param.file '\'];
    Zim = imreadstack([directory param.Date '_' param.file '_MMStack_Default.ome.tif']);
    im = imreadstack([directory param.Date '_' param.file '_MMStack_Default_561_denoise.ome.tif']);
    tracks = xlsread([directory param.Date '_' param.file '_tracks.csv']);
    cells = xlsread([directory param.Date '_' param.file '_MMStack_Default_561_denoise.ome.roi.zip.csv']);
    
    param.derivative_threshold = 400;
    param.pixSz = 65; % [nm/pix]
    param.interval = 1; % [s]
    param.psfFWHM = 250; % [nm]
else
    if isempty(directory)
        path = uigetdir;
    else
        path = directory;
    end
    if nargin < 4
        param.derivative_threshold = 40;
    end
    if nargin < 3
        param.pixSz = 65; % [nm/pix]
    end
    
    Zim_file = uigetfile('*.ome.tif', 'FtsZ image file', [path '\']);
    im_file = uigetfile('*.ome.tif', 'Pbp2B stack', [path '\']);
    tracks_file = uigetfile('*.csv', 'Track file', [path '\']);
    septa_file = uigetfile('*.ome.roi.zip.csv', 'Septa file', [path '\']);
    
    Zim = imreadstack([path '\' Zim_file]);
    im = imreadstack([path '\' im_file]);
    tracks = xlsread([path '\' tracks_file]);
    cells = xlsread([path '\' septa_file]);
end

% USER INPUTS

plot_fov = 0; % plot the field of view with tracks and septa
plot_vel_msd = 0; % plot a histogram of velocities obtaines from MSD analysis

pix2um = param.pixSz/1000; % [um/pix]

% filters
param.time_range = [10 1e10] / param.interval; % [frames] threshold for how long a track needs to last
param.pairing_distance_threshold = 0.5; % [um] threshold for how far a track center can be from a septum center
param.distance_threshold = 0.05; % [um] threshold for end-to-end distance a track needs to go
param.R2_threshold = 0.95; % threshold for R^2 (for fit of log(MSD) to log(t))

% no filters
param.distance_threshold = 0.00; % [um] threshold for end-to-end distance a track needs to go
param.R2_threshold = 0.0; % threshold for R^2 (for fit of log(MSD) to log(t))

% sort rows (times sometimes are scrambled)
tracks = sortrows(tracks,[2 7]);

% Start plots

if plot_fov
    figure
    ax_xy = gca;
    hold on
    ylabel('y (\mum)')
    xlabel('x (\mum)')
    set(ax_xy,"YDir","reverse")
    xlim([0 66.56])
    ylim([0 66.56])
    axis equal
end

if plot_vel_msd
    figure
    ax_vel_msd = gca;
    hold on
    ylabel('Counts')
    xlabel('Speed from MSD fit (nm/s)')
end

% make cell array of septal (i.e. cell) coordinates

% septa(:,3:4) = septa(:,3:4).*pix2um; % [um] septal coordinates
cells = sortrows(cells, [1 2]); % rearrange so first x coordinate is always lowest, second is highest

Cells = {}; cc1 = 1;
for cc = 1:4:length(cells) % data are by fours, so increment by 4

    Cells{cc1}(:,1) = [cells(cc,3); cells(cc+1,3)];
    Cells{cc1}(:,2) = [cells(cc+1,4); cells(cc+2,4)];
    
%     plot(p_xy, Septa{ss1}(:,1), Septa{ss1}(:,2))
    
    cc1 = cc1+1;
end

% separate tracks

tracknums = unique(tracks(:,2));

Velocities=[]; dr=[]; subim=zeros(60,60);
for ii = 1:length(Cells)
    
    ud.cell.cellnum = ii;
    
    % start plot for this cell
    ud.figsavename = [directory param.Date '_' param.file '_cell' num2str(ii)];
    fig_tracks = figure('Position',[850, 130, 770, 800], 'FileName', ud.figsavename);
    ud.axes.ax_im = subplot(3,4,1);
    hold on
    axis equal
    
    ud.axes.ax_polar = subplot(3,3,2,polaraxes);
    hold on
    filename_plot = replace([param.Date '_' param.file], '_', '\_');
    title(ud.axes.ax_polar, ['HaloTag-PBP2B tracking data for FoV ' filename_plot ', cell ' num2str(ii)])
    ud.axes.ax_polar.ThetaDir = 'clockwise';
    ud.axes.ax_polar.ThetaTick = 0:90:360;
    
    ud.axes.ax_msd = subplot(3,3,3);
    hold on
    ylabel('MSD (\mum^2)')
    xlabel('\Deltat (s)')
    set(ud.axes.ax_msd,'YAxisLocation','right')
    
    ud.axes.ax_rho = subplot(3,3,4:6);
    hold on
    ylabel('\rho (nm)')
    ud.axes.ax_dist = subplot(3,3,7:9);
    hold on
    ylabel('Arc length (nm)')
    xlabel('Time (s)')
    
    % FIND COORDINATES FOR CENTER OF RING USING Z-RING IMAGE
    
    % create bounding box (need to cut if near border)
    box_x1 = max([Cells{ii}(1,1), 1]); % [pix]
    box_x2 = min([Cells{ii}(2,1), size(im,1)]); % [pix]
    box_y1 = max([Cells{ii}(1,2), 1]); % [pix]
    box_y2 = min([Cells{ii}(2,2), size(im,2)]); % [pix]
    xdim = box_x2 - box_x1; % [pix]
    ydim = box_y2 - box_y1; % [pix]
    
    ud.cell.boxsize = [xdim ydim];
    
    % get individual ring bounded by the box at index ind_box
%     subim(1:ydim+1,1:xdim+1,1:size(im,3)) = im(box_y1:box_y2, box_x1:box_x2, :); % need to flip x and y
    subim(1:ydim+1,1:xdim+1,1) = Zim(box_y1:box_y2, box_x1:box_x2, 1);
    
    % fit individual ring to model, get center position
    %ringIm = mean(subim,3);
%     ringImMax = max(subim,[],3);
    fitRingArg = {'FixedPositionFit', true};
    [fitParAvg, ~, im_bg_sub] = fitRing_public(subim, pix2um*1000, param.psfFWHM, fitRingArg{:});
    cent = fitParAvg(:,1:2); % [pix]
    imagesc(ud.axes.ax_im, im_bg_sub)
    plot(ud.axes.ax_im, cent(1), cent(2), 'xr')
    tv_theta = 0:0.01:2*pi;
    plot(ud.axes.ax_im, cent(1)-fitParAvg(:,3)*cos(tv_theta), cent(2)-fitParAvg(:,3)*sin(tv_theta),'r')
    
    % plot square (user-chosen in ImageJ) and septum center from fits
    if plot_fov
        rectangle(ax_xy, 'Position', [Cells{ii}(1,1) Cells{ii}(1,2) 60 60].*pix2um)
        plot(ax_xy, (cent(1)+Cells{ii}(1,1))*pix2um, (cent(2)+Cells{ii}(1,2))*pix2um, 'xk')
    end
    
    ud.cell.cell_center_fit = cent;
    ud.cell.cell_radius_fit = fitParAvg(:,3);
    
    % PAIR SEPTUM WITH TRACK
    
    % pair septum with track
    mean_septa = mean(Cells{ii},1).*pix2um;
    for ss2 = 1:length(tracknums)
        median_track_pos = median(tracks(tracks(:,2)==tracknums(ss2),4:5));
        dr(ss2) = sqrt((median_track_pos(1)-mean_septa(1)).^2 + (median_track_pos(2)-mean_septa(2)).^2); % distance difference
    end
    
    % find tracks associated with that cell
    tracknum_paired = tracknums(dr <= param.pairing_distance_threshold);
    
    ud.cell.Ntracks = length(tracknum_paired);
    ud.cell.tracknums = tracknum_paired;

    % go through each track for this cell
    Tracks={}; maxrho=[]; ud.phandles.p_polar=[]; ud.phandles.p_rho=[]; ud.phandles.p_dist=[]; ud.phandles.p_msd=[]; ud.phandles.p_msd_fit=[];
    for tt = 1:length(tracknum_paired)
        
        ud.cell.track{tt}.trackind = tt;
        ud.cell.track{tt}.tracknum = tracknum_paired(tt);
        
        Tracks{tt} = tracks(tracks(:,2)==tracknum_paired(tt),:);
        
        % plot track with associated cell
        if plot_fov
            plot(ax_xy, Tracks{tt}(:,4), Tracks{tt}(:,5)) % track
        end
        
        % REMOVE TRACKS THAT ARE TOO SHORT OR TOO LONG
        
        if size(Tracks{tt},1) < param.time_range(1) || size(Tracks{tt},1) > param.time_range(2)
            continue
        end
        
        
        % REMOVE TRACKS THAT DIDN'T GO FAR
        
        Rxy = Tracks{tt}(:,4:5); % [um] (x,y) pairs
        
        if norm(Rxy(end,:)-Rxy(1,:)) < param.distance_threshold
            continue
        end
        
        
        % CONVERT TRACK DATA TO POLAR COORDINATES
        
        % subtract center position from the track
        Tracks{tt}(:,20) = (cent(1)+box_x1)*pix2um - Tracks{tt}(:,4); % [um]
        Tracks{tt}(:,21) = (cent(2)+box_y1)*pix2um - Tracks{tt}(:,5); % [um]
        Tracks{tt}(:,20) = -Tracks{tt}(:,20); % reflect about x
        Tracks{tt}(:,21) = -Tracks{tt}(:,21); % reflect about y. needed because we want to clockwise polar plot.
        
        % get polar coordinates
        Tracks{tt}(:,22) = sqrt(Tracks{tt}(:,20).^2 + Tracks{tt}(:,21).^2) * 1000; % [nm] rho
        Tracks{tt}(:,23) = atan2(Tracks{tt}(:,21), Tracks{tt}(:,20)); % [rad] theta
        
        ud.cell.track{tt}.time = Tracks{tt}(:,7)*param.interval;
        ud.cell.track{tt}.rho = Tracks{tt}(:,22);
        ud.cell.track{tt}.theta = Tracks{tt}(:,23);
        ud.cell.track{tt}.rho_new = Tracks{tt}(:,22); % duplicate in case you shift it later
        ud.cell.track{tt}.theta_new = Tracks{tt}(:,23); % duplicate in case you shift it later
        
        Tracks{tt}(:,24) = Tracks{tt}(:,23).*Tracks{tt}(:,22); % [nm] angular position
        Tracks{tt}(:,23) = mod(Tracks{tt}(:,23), 2*pi); % [rad] convert from [-pi, pi] to [0, 2*pi]
        Tracks{tt}(:,23) = Tracks{tt}(:,23) * 180/pi; % [deg]
        
        % patch up boundary conditions
        dxdt = (Tracks{tt}(2:end,24)-Tracks{tt}(1:end-1,24)) / param.interval; % first derivative
        ind_plus = find(dxdt>=param.derivative_threshold);
        ind_minus = find(dxdt<=-param.derivative_threshold);
        
        % when track reaches boundary [-pi pi], add or subtract 2*pi*r to keep track continuous
        for jj = 1:length(ind_plus)
            Tracks{tt}(ind_plus(jj)+1:end,24) = Tracks{tt}(ind_plus(jj)+1:end,24) - 2*pi*Tracks{tt}(ind_plus(jj)+1:end,22);
        end
        for jj = 1:length(ind_minus)
            Tracks{tt}(ind_minus(jj)+1:end,24) = Tracks{tt}(ind_minus(jj)+1:end,24) + 2*pi*Tracks{tt}(ind_minus(jj)+1:end,22);
        end
        
        % REMOVE TRACKS WITH POOR FITS TO MSD VS. TIME
        
        [RX1, RX2] = meshgrid(Rxy(:,1),Rxy(:,1));
        [RY1, RY2] = meshgrid(Rxy(:,2),Rxy(:,2));
        
        SD = (RX2-RX1).^2 + (RY2-RY1).^2; % [um^2] squared displacements
        SD = triu(SD); % [um^2] upper triangular part of matrix
        SD = spdiags(SD); % [um^2] rearrange off-diagonals into columns
        SD(SD==0) = nan; % [um^2] remove zeros
        MSD = nanmean(SD,1); % [um^2] mean squared displacements
        
        Dt = param.interval:param.interval:length(MSD)*param.interval; % [s] time delays
        n_MSD = sum(~isnan(SD)); % number of MSDs for each time delay
        
        eqn = @(a,x) a(1) + a(2).*x; % a(1) = log(4D), a(2) = alpha, x = log(t)
        [fitvals, ~, ~] = nlinfit(log(Dt), log(MSD), eqn, [1e-3 1], 'Weights', n_MSD);
        
        sse = sum((log(MSD)-eqn(fitvals,log(Dt))).^2); % residual sum of squares
        sst = sum((log(MSD)-mean(log(MSD))).^2); % total sum of squares
        r2 = 1 - sse/sst; % R squared
        
        if r2 < param.R2_threshold
            continue
        end
        
        % CALCULATE VELOCITY FROM FITS TO MSD
        
        eqn_vel = @(a,x) 4*a(1)*x + (a(2)*x).^2; % a(1) = D, a(2) = V
        fit_vel = nlinfit(Dt, MSD, eqn_vel, [1e-6 0.01], 'Weights', n_MSD);
        
        ud.cell.track{tt}.msd_eqn = eqn_vel;
        ud.cell.track{tt}.msd_D = fit_vel(1);
        ud.cell.track{tt}.msd_V = fit_vel(2);
        
        Velocities = [Velocities; fit_vel(2)];

        
        % PLOT TRACKS, MSD, AND TRAJECTORIES
        
        % make polar plot
        ud.phandles.p_polar(tt) = polarplot(ud.axes.ax_polar, Tracks{tt}(:,23)*pi/180, Tracks{tt}(:,22), 'DisplayName', ['Track ' num2str(tt)]);
        track_clr = get(ud.phandles.p_polar(tt),'color');
        
        % plot MSD
        %     plot(p_msd, log(Dt), log(MSD), 'ok')
        %     plot(p_msd, log(Dt), eqn(fitvals,log(Dt)), 'r')
        ud.phandles.p_msd(tt) = plot(ud.axes.ax_msd, Dt, MSD, '.', 'Color', track_clr, 'DisplayName', ['Track ' num2str(tt)]);
        ud.phandles.p_msd_fit(tt) = plot(ud.axes.ax_msd, Dt, eqn_vel(fit_vel,Dt), 'Color', track_clr, 'DisplayName', ['Track ' num2str(tt)]);
        %     plot(p_msd, Dt, MSDpll, 'ob')
        %     plot(p_msd, Dt, eqn1D(fitpll,Dt), 'g')
        
        % plot in polar coordinates
        ud.phandles.p_rho(tt) = plot(ud.axes.ax_rho, Tracks{tt}(:,7)*param.interval, Tracks{tt}(:,22), 'DisplayName', ['Track ' num2str(tt)]); % rho
        maxrho(tt) = max(Tracks{tt}(:,22));
        ud.phandles.p_dist(tt) = plot(ud.axes.ax_dist, Tracks{tt}(:,7)*param.interval, Tracks{tt}(:,24)-Tracks{tt}(1,24), 'DisplayName', ['Track ' num2str(tt)]); % theta
    end
    
    if isempty(ud.axes.ax_polar.Children) %  no tracks plotted - delete figure 
        delete(fig_tracks)
        continue
    end
    
    ud.axes.ax_msd.XScale = 'log';
    ud.axes.ax_msd.YScale = 'log';
    ud.axes.ax_rho.YLim = [0 max(maxrho)*1.1];
    linkaxes([ud.axes.ax_rho ud.axes.ax_dist],'x')
    
    ud.param = param;
    fig_tracks.UserData = ud;
    
    set(fig_tracks,'WindowKeyPressFcn',@multikey_analysis)

end

if plot_vel_msd
    histogram(ax_vel_msd, Velocities.*1000, 0:2:60)
end
