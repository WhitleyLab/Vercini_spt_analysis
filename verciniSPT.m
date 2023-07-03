
function verciniSPT(param)

if nargin == 0
    param.Date = '220316';
    param.file = '1';
    param.analysisDate = datestr(now,'yymmdd');
    param.directory = ['\\campus\rdw\FMS CBCB\nsh167\Shared\data\Whitley-Kevin\220316_sh147_phmm_30c_vercini\' param.Date '_' param.file '\'];
%     Zim = imreadstack([param.directory param.Date '_' param.file '_MMStack_Default.ome.tif']);
    im = imreadstack([param.directory param.Date '_' param.file '_MMStack_Default_561_denoise.ome.tif']);
    tracks = xlsread([param.directory param.Date '_' param.file '_tracks_diffusive_nogaps.csv']);
    cells = xlsread([param.directory param.Date '_' param.file '_MMStack_Default_561_denoise.ome.roi.zip.csv']);

    param.pixSz = 65; % [nm/pix]
    param.interval = 1; % [s]
    param.psfFWHM = 250; % [nm]
    param.derivative_threshold = 1000 / param.interval;
    param.fixedRadius = 1; % use fixed radius (otherwise use instantaneous 'radius')
    param.thunderstorm = 0; % used thunderSTORM for localizations and tracking (tracking csv file in um for regular TrackMate, nm for tSTORM+TrackMate)
    
    param.diffusion_model = 'brownian'; % brownian or directional
    param.msd_fit = 'loglog'; % use loglog fit or linear fit?
    param.msd_fit_range = [0 15]; % range of time-steps to fit for MSD data
else
    if isempty(param.directory)
        path = uigetdir;
    else
        path = param.directory;
    end
    if nargin < 4
        param.derivative_threshold = 40;
    end
    if nargin < 3
        param.pixSz = 65; % [nm/pix]
    end
    
%     Zim_file = uigetfile('*.ome.tif', 'FtsZ image file', [path '\']);
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
param.time_range = [10 Inf]; % [frames] threshold for how long a track needs to last
param.pairing_distance_threshold = 2; % [um] threshold for how far a track center can be from a septum center
param.distance_threshold = 0.05; % [um] threshold for end-to-end distance a track needs to go
param.R2_threshold = 0.95; % threshold for R^2 (for fit of log(MSD) to log(t))

% no filters
param.distance_threshold = 0.00; % [um] threshold for end-to-end distance a track needs to go
param.R2_threshold = -Inf; % threshold for R^2 (for fit of log(MSD) to log(t))

% if (x,y) coordinates in track file were in nm, convert to um
if param.thunderstorm
    tracks(:,4:5) = tracks(:,4:5)./1000; % [um]
    tracks(:,7) = tracks(:,8); % [frames] columns mixed up in this format
end

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

Velocities=[]; dr=[]; subim=zeros(60,60); subim_2b=zeros(60,60,size(im,3));
for ii = 1:length(Cells)
    
    ud.cell.cellnum = ii;
    
    % start plot for this cell
    ud.figsavename = [param.directory param.Date '_' param.file '_cell' num2str(ii) '_diffuse'];
    fig_tracks = figure('Position',[850, 130, 770, 800], 'FileName', ud.figsavename);
%     ud.axes.ax_im = subplot(3,4,1);
    ud.axes.ax_im = subplot(4,4,1);
    set(ud.axes.ax_im,"YDir","reverse")
    hold on
    axis equal
    
%     ud.axes.ax_polar = subplot(3,3,2,polaraxes);
    ud.axes.ax_polar = subplot(4,3,2,polaraxes);
    hold on
    filename_plot = replace([param.Date '_' param.file], '_', '\_');
    title(ud.axes.ax_polar, ['HaloTag-PBP2B tracking data for FoV ' filename_plot ', cell ' num2str(ii)])
    ud.axes.ax_polar.ThetaDir = 'clockwise';
    ud.axes.ax_polar.ThetaTick = 0:90:360;
    
%     ud.axes.ax_msd = subplot(3,3,3);
    ud.axes.ax_msd = subplot(4,3,3);
    hold on
    ylabel('MSD (\mum^2)')
    xlabel('\Deltat (s)')
    set(ud.axes.ax_msd,'YAxisLocation','right')
    
%     ud.axes.ax_rho = subplot(3,3,4:6);
    ud.axes.ax_rho = subplot(4,3,4:6);
    hold on
    ylabel('\rho (nm)')
%     ud.axes.ax_dist = subplot(3,3,7:9);
    ud.axes.ax_dist = subplot(4,3,7:9);
    hold on
    ylabel('Arc length (nm)')
    
    ud.axes.ax_int = subplot(4,3,10:12);
    hold on
    ylabel('Intensity')
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
    subim(1:ydim+1,1:xdim+1,1:size(im,3)) = im(box_y1:box_y2, box_x1:box_x2, :); % need to flip x and y
%     subim(1:ydim+1,1:xdim+1,1) = Zim(box_y1:box_y2, box_x1:box_x2, 1);
    
    % fit individual ring to model, get center position
    %ringIm = mean(subim,3);
    ringImMax = max(subim,[],3);
    fitRingArg = {'FixedPositionFit', true};
    try
        [fitParAvg, ~, im_bg_sub] = fitRing_public(ringImMax, pix2um*1000, param.psfFWHM, fitRingArg{:});
    catch
        continue
    end
%     [fitParAvg, ~, im_bg_sub] = fitRing_public(subim, pix2um*1000, param.psfFWHM, fitRingArg{:});
    cent = fitParAvg(:,1:2); % [pix]
    radius = fitParAvg(:,3); % [pix]
    implot = imagesc(ud.axes.ax_im, im_bg_sub);
    set(implot, 'XData', implot.XData-0.5) % make sure pixels actually start at (0,0)
    set(implot, 'YData', implot.YData-0.5) % make sure pixels actually start at (0,0)
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
    
    
    % GET BACKGROUND FLUORESCENCE
    
    subim_2b(1:ydim+1,1:xdim+1,:) = im(box_y1:box_y2, box_x1:box_x2, :);
    imBlur = imgaussfilt(subim_2b,1);
    bgval=[]; otsuThresh=[];
    for kk = 1:size(im,3)
       otsuThresh(:,:,kk) = otsu_public(imBlur(:,:,ii));
        oneframe = subim_2b(:,:,kk);
        bgval(kk) = median(oneframe(otsuThresh(:,:,kk)==1));
    end

    eqn_exp = @(a,x) a(1)*exp(-x/a(2)) + a(3);
    cell_time = 0:param.interval:(size(subim_2b,3)-1)*param.interval;
    
    init_bg_fit = [max(bgval)-min(bgval) size(subim_2b,3)/2 min(bgval)];
    intFit = nlinfit(cell_time, bgval, eqn_exp, init_bg_fit);
    
    bg = eqn_exp(intFit, cell_time); % exponential decay in background signal over time

    
    % GO THROUGH EACH TRACK

    % go through each track for this cell
    Tracks={}; maxrho=[]; ud.phandles.p_polar=[]; ud.phandles.p_rho=[]; ud.phandles.p_dist=[]; ud.phandles.p_msd=[]; ud.phandles.p_msd_fit=[]; ud.cell.track={}; JD_all=[];
    for tt = 1:length(tracknum_paired)
        
        ud.cell.track{tt}.trackind = tt;
        ud.cell.track{tt}.tracknum = tracknum_paired(tt);
        ud.cell.track{tt}.isplotted = 0; % changes to 1 if this track is plotted
        
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
        cellcent_fov_x = (cent(1)+box_x1-1)*pix2um; % [um] x position of center of cell relative to full FoV. included -1 (220613 kw)
        cellcent_fov_y = (cent(2)+box_y1-1)*pix2um; % [um] y position of center of cell relative to full FoV. included -1 (220613 kw)
        Tracks{tt}(:,20) = cellcent_fov_x - Tracks{tt}(:,4); % [um]
        Tracks{tt}(:,21) = cellcent_fov_y - Tracks{tt}(:,5); % [um]
        Tracks{tt}(:,20) = -Tracks{tt}(:,20); % reflect about x
        Tracks{tt}(:,21) = -Tracks{tt}(:,21); % reflect about y. needed because we want to clockwise polar plot.
        
        % get polar coordinates
        Tracks{tt}(:,22) = sqrt(Tracks{tt}(:,20).^2 + Tracks{tt}(:,21).^2) * 1000; % [nm] rho
        Tracks{tt}(:,23) = atan2(Tracks{tt}(:,21), Tracks{tt}(:,20)); % [rad] theta
        
        ud.cell.track{tt}.time = Tracks{tt}(:,7) * param.interval; % [s]
        if param.thunderstorm % HACK 220622. don't know why I have to do this. format is different again.
            ud.cell.track{tt}.intensity = Tracks{tt}(:,12) - bg(Tracks{tt}(:,7))';
        else
            ud.cell.track{tt}.intensity = Tracks{tt}(:,13) - bg(Tracks{tt}(:,7)+1)';
        end
        
        ud.cell.track{tt}.rho = Tracks{tt}(:,22);
        ud.cell.track{tt}.theta = Tracks{tt}(:,23);
        ud.cell.track{tt}.rho_new = Tracks{tt}(:,22); % duplicate in case you shift it later
        ud.cell.track{tt}.theta_new = Tracks{tt}(:,23); % duplicate in case you shift it later
        
        if param.fixedRadius
            Tracks{tt}(:,24) = Tracks{tt}(:,23) .* radius*pix2um*1000; % [nm] angular position, fixed R
        else
            Tracks{tt}(:,24) = Tracks{tt}(:,23) .* Tracks{tt}(:,22); % [nm] angular position, instantaneous R
        end
        Tracks{tt}(:,23) = mod(Tracks{tt}(:,23), 2*pi); % [rad] convert from [-pi, pi] to [0, 2*pi]
        Tracks{tt}(:,23) = Tracks{tt}(:,23) * 180/pi; % [deg]
        
        % patch up boundary conditions
        dxdt = (Tracks{tt}(2:end,24)-Tracks{tt}(1:end-1,24)) / param.interval; % first derivative
        ind_plus = find(dxdt>=param.derivative_threshold);
        ind_minus = find(dxdt<=-param.derivative_threshold);
        
        % when track reaches boundary [-pi pi], add or subtract 2*pi*r to keep track continuous
        for jj = 1:length(ind_plus)
            if param.fixedRadius
                Tracks{tt}(ind_plus(jj)+1:end,24) = Tracks{tt}(ind_plus(jj)+1:end,24) - 2*pi*radius*pix2um*1000;
            else
                Tracks{tt}(ind_plus(jj)+1:end,24) = Tracks{tt}(ind_plus(jj)+1:end,24) - 2*pi*Tracks{tt}(ind_plus(jj)+1:end,22);
            end
        end
        for jj = 1:length(ind_minus)
            if param.fixedRadius
                Tracks{tt}(ind_minus(jj)+1:end,24) = Tracks{tt}(ind_minus(jj)+1:end,24) + 2*pi*radius*pix2um*1000;
            else
                Tracks{tt}(ind_minus(jj)+1:end,24) = Tracks{tt}(ind_minus(jj)+1:end,24) + 2*pi*Tracks{tt}(ind_minus(jj)+1:end,22);
            end
        end
        
        
        % CALCULATE JUMP DISTANCES FOR EACH Dt (both 1D and 2D)
        
        % 1D case
        
        [RS1, RS2] = meshgrid(Tracks{tt}(:,24), Tracks{tt}(:,24)); % [nm] grid of arc lengths
        
        JD_1D = abs(RS2-RS1) ./ 1000; % [um] jump distances along circumference (1D)
        JD_1D = triu(JD_1D); % [um] upper triangular part of matrix
        JD_1D = spdiags(JD_1D); % [um] rearrange off-diagonals into columns (each column now a different Dt)
        JD_1D(JD_1D==0) = nan; % [um] remove zeros
        
        ud.cell.track{tt}.jump_distances_1D = JD_1D;
        
        % 2D case
        
        [RX1, RX2] = meshgrid(Rxy(:,1),Rxy(:,1));
        [RY1, RY2] = meshgrid(Rxy(:,2),Rxy(:,2));
        
        JD_2D = sqrt((RX2-RX1).^2 + (RY2-RY1).^2); % [um] jump distances
        JD_2D = triu(JD_2D); % [um] upper triangular part of matrix
        JD_2D = spdiags(JD_2D); % [um] rearrange off-diagonals into columns (each column now a different Dt)
        JD_2D(JD_2D==0) = nan; % [um] remove zeros
        
        ud.cell.track{tt}.jump_distances_2D = JD_2D;

        
        % CALCULATE MEAN SQUARED DISPLACEMENTS AND FIT TO DIFFUSION MODEL

        SD = (RX2-RX1).^2 + (RY2-RY1).^2; % [um^2] squared displacements
        SD = triu(SD); % [um^2] upper triangular part of matrix
        SD = spdiags(SD); % [um^2] rearrange off-diagonals into columns
        SD(SD==0) = nan; % [um^2] remove zeros
        MSD = nanmean(SD,1); % [um^2] mean squared displacements
        
        Dt = param.interval:param.interval:length(MSD)*param.interval; % [s] time delays
        n_MSD = sum(~isnan(SD)); % number of MSDs for each time delay
        
        fit_ind = min([param.msd_fit_range(2) length(Dt)]);
        Dt_cut = Dt(1:fit_ind);
        MSD_cut = MSD(1:fit_ind);
        n_MSD_fit = n_MSD(1:fit_ind);
        
        % Choose diffusion model
        if strcmp(param.diffusion_model, 'brownian') % anomalous diffusion model (fit to D and alpha)
            if strcmp(param.msd_fit, 'loglog')
                eqn_diff = @(a,x) a(1) + x.*a(2); % log version. a(1) = log(2n*D), a(2) = alpha, where n is dimensionality
                Dt_fit = log(Dt_cut); % for fitting, use log of Dt
                MSD_fit = log(MSD_cut); % for fitting, use log of MSD
                init_diff = [-11 1];
                Dt_plot = log(Dt); % for plotting the model, use log of Dt
            else
                eqn_diff = @(a,x) a(1)*x.^a(2); % a(1) = 2n*D, a(2) = alpha, where n is dimensionality
                Dt_fit = Dt_cut; % for fitting, use log of Dt
                MSD_fit = MSD_cut; % for fitting, use log of MSD
                init_diff = [1e-3 1];
                Dt_plot = Dt; % for plotting the model
            end
                %             eqn_diff = @(a,x) a(1) + log(x).*a(2); % log version. a(1) = log(2n*D), a(2) = alpha, where n is dimensionality
            
%             eqn_fit = fittype(@(a,b,x)eqn_diff(a,b,x));
            
%             init_diff = [-11 1];
%             lb = [-Inf 0];
%             ub = [Inf Inf];
%             options = optimoptions('lsqcurvefit', 'Display', 'off');
%             [fitvals, ~, res, ~, ~, ~, J] = lsqcurvefit(eqn_diff, init_diff, Dt_fit, MSD_fit, lb, ub, options);
            [fitvals, res, J] = nlinfit(Dt_fit, MSD_fit, eqn_diff, init_diff);
            ci = nlparci(fitvals, res, 'jacobian', J);
            t = tinv(1-0.05/2, length(MSD_fit)-length(fitvals));
            fiterrs = (ci(:,2)-ci(:,1)) ./ (2*t);
        else % diffusion + directional motion model (fit to D and v)
            if strcmp(param.msd_fit, 'loglog') % actually in this case it's semilog
                eqn_diff = @(a,x) log(a(1)*x + (a(2)*x).^2); % log version. a(1) = 2n*D, a(2) = V, where n is dimensionality
                MSD_fit = log(MSD_cut); % for fitting, use log of MSD
                init_diff = [0.05 0.01];
            else % linear fit
                eqn_diff = @(a,x) a(1)*x + a(2)*x.^2; % a(1) = 2n*D, a(2) = V, where n is dimensionality
                MSD_fit = MSD_cut;
                init_diff = [1e-3 0.01];
            end
%             eqn_fit = fittype(@(a,b,x)eqn_diff(a,b,x));
            
            Dt_fit = Dt_cut; % for fitting, use regular Dt (not log)

%             lb = [0 0];
%             ub = [Inf Inf];
%             options = optimoptions('lsqcurvefit', 'Display', 'off');
%             fitvals = lsqcurvefit(eqn_diff, init_diff, Dt_fit, MSD_fit, lb, ub, options);
            [fitvals, res, J] = nlinfit(Dt_fit, MSD_fit, eqn_diff, init_diff);
            ci = nlparci(fitvals, res, 'jacobian', J);
            t = tinv(1-0.05/2, length(MSD_fit)-length(fitvals));
            fiterrs = (ci(:,2)-ci(:,1)) ./ (2*t);
            
            Dt_plot = Dt; % fot plotting the model, use regular Dt (not log)
        end

%         [fitvals, gof] = fit(Dt_fit', MSD_fit', eqn_diff, 'Weights', n_MSD_fit, 'StartPoint', init_diff);
%         [fitvals, ~, ~] = nlinfit(Dt_fit, log(MSD_fit), eqn_diff, init_diff, 'Weights', n_MSD_fit);

        sse = sum((log(MSD)-eqn_diff(fitvals,log(Dt))).^2); % residual sum of squares
        sst = sum((log(MSD)-mean(log(MSD))).^2); % total sum of squares
        r2 = 1 - sse/sst; % R squared
        
        if r2 < param.R2_threshold % filter by goodness of fit to MSD
            continue
        end
        
        ud.cell.track{tt}.diffusion_eqn = eqn_diff;
        ud.cell.track{tt}.msd_fits = fitvals; % [D alpha] or [D v], depending on which model you chose
        ud.cell.track{tt}.msd_errs = fiterrs; % [D alpha] or [D v], depending on which model you choose
        
        
        % PLOT TRACKS, MSD, AND TRAJECTORIES
        
        % make polar plot
        ud.phandles.p_polar(tt) = polarplot(ud.axes.ax_polar, Tracks{tt}(:,23)*pi/180, Tracks{tt}(:,22), 'DisplayName', ['Track ' num2str(tt)]);
        track_clr = get(ud.phandles.p_polar(tt),'color');
        ud.cell.track{tt}.isplotted = 1;
        
        % plot MSD
        ud.phandles.p_msd(tt) = plot(ud.axes.ax_msd, Dt, MSD, '.', 'Color', track_clr, 'DisplayName', ['Track ' num2str(tt)]);
        if strcmp(param.msd_fit, 'loglog')
            ud.phandles.p_msd(tt) = plot(ud.axes.ax_msd, Dt, exp(eqn_diff(fitvals,Dt_plot)), 'Color', track_clr, 'DisplayName', ['Track ' num2str(tt)]);
        else
            ud.phandles.p_msd(tt) = plot(ud.axes.ax_msd, Dt, eqn_diff(fitvals,Dt_plot), 'Color', track_clr, 'DisplayName', ['Track ' num2str(tt)]);
        end
%         ud.phandles.p_msd(tt) = plot(ud.axes.ax_msd, Dt, eqn_diff(fitvals.a, fitvals.b, Dt), 'Color', track_clr, 'DisplayName', ['Track ' num2str(tt)]);
        
        % plot in polar coordinates
        ud.phandles.p_rho(tt) = plot(ud.axes.ax_rho, Tracks{tt}(:,7)*param.interval, Tracks{tt}(:,22), 'DisplayName', ['Track ' num2str(tt)]); % rho
        maxrho(tt) = max(Tracks{tt}(:,22));
%         ud.phandles.p_dist(tt) = plot(ud.axes.ax_dist, Tracks{tt}(:,7)*param.interval, Tracks{tt}(:,24)-Tracks{tt}(1,24), 'DisplayName', ['Track ' num2str(tt)]); % theta
        ud.phandles.p_dist(tt) = plot(ud.axes.ax_dist, Tracks{tt}(:,7)*param.interval, Tracks{tt}(:,24), 'DisplayName', ['Track ' num2str(tt)]); % theta
        
        ud.phandles.p_int(tt) = plot(ud.axes.ax_int, Tracks{tt}(:,7)*param.interval, ud.cell.track{tt}.intensity, 'DisplayName', ['Track ' num2str(tt)]); % intensity
        maxint(tt) = max(ud.cell.track{tt}.intensity);
    end
    
    if isempty(ud.axes.ax_polar.Children) %  no tracks plotted - delete figure 
        delete(fig_tracks)
        continue
    end
    
    ud.axes.ax_msd.XScale = 'log';
    ud.axes.ax_msd.YScale = 'log';
    ud.axes.ax_rho.YLim = [0 max(maxrho)*1.1];
    ud.axes.ax_dist.YLim = [-pi*radius*pix2um*1000 pi*radius*pix2um*1000];
    ud.axes.ax_dist.XLim = [0 10];
    ud.axes.ax_int.YLim = [0 max(maxint)*1.1];
    linkaxes([ud.axes.ax_rho ud.axes.ax_dist ud.axes.ax_int],'x')
    
    ud.param = param;
    fig_tracks.UserData = ud;
    
    set(fig_tracks,'WindowKeyPressFcn',@multikey_analysis)

end

if plot_vel_msd
    histogram(ax_vel_msd, Velocities.*1000, 0:2:60)
end
