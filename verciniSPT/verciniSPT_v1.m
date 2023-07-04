% Author: Kevin Whitley
% 
% This function plots single-particle tracking data for vertically-trapped
% bacterial cells. It takes a video file (.tif format), an associated
% tracking data file (.csv format), and an associated data file with the
% coordinates of cells to be analyzed (.csv format). It then plots the
% tracks in polar coordinates, measures and fits the MSDs, and plots the
% time traces in polar coordinates. The figure is then interactive, with
% several options for further analysis.
% 
% FILE FORMATS:
%   Video file: Should be in .tif format
%   Track data file: Contains output of TrackMate data, in .csv format
%   Cell coordinate file: Contains coordinates of the bounding boxes of cells chosen manually in ImageJ, in .csv format
% 
% INPUTS:
%   fileFilter: Search string for files to analyze (e.g. '*.tif')
%   interval: Acquisition interval [s]
%   pixSz: Camera pixel size [nm/px]
%   varargin: Optional inputs:
%       'PlotFOV', true/false (default:false): Plot the Field of View with all tracks and bounding boxes.
%       'PsfFWHM', (default:250): FWHM of the point-spred function [nm].
%       'TimeRange', (default:[10 Inf]): Time range to use for filtering tracks (only tracks in this range will be used for subsequent analysis [frames].
%       'DistThresh', (default:0): Threshold for how far a track needs to go (end-to-end) to be included [um]. Only tracks with end-to-end distance greater than this will be used in subsequent analysis.
%       'PairDistThresh', (default:2): Threshold for how far a track center can be from a septum center [um]. Only tracks within this distance from the calculated center of a track will be used in subsequent analysis.
%       'DiffusionModel', ('brownian', 'directional') (default:'brownian'): Which model to use for fits to MSD vs. timestep plots. 'brownian' uses the anomalous diffusion model MSD = 4Dt^alpha, while 'directional' uses MSD = 4Dt + (vt)^2.
%       'MSDfit', ('linear', 'loglog') (default:'loglog'): Fit MSD vs. timestep plots on linear scale or on log-log scale.
%       'MSDfitrange', (default:[0 15]): Time range to use for fitting MSD vs. timestep plots [s].
%       'RsquaredThresh', (default:-Inf): Threshold for R^2 for MSD vs. timestep plots. Only tracks with R^2 greater than this will be used in subsequent analysis.
% 
% OUTPUTS:
%   A figure file with:
%       An image of the cell being analyzed with fits to a ring.
%       A plot of all tracks analyzed in polar coordinates.
%       A plot of MSDs vs. timestep with fits for each track.
%       A time trace of coordinate rho (distance from center) for each track.
%       A time trace of arc length (position around circumference) for each track.
%       A time trace of intensity for each track.
%   Figures have saved within them all parameters used for analysis and results for any fits (e.g. MSD fits) as the figure's UserData.
% 
%   Figures are interactive. They can be further analyzed by the user by pressing specific keys and clicking points on the figure:
%       Press 'c' and click on polar plot to reset the center position of the cell. (often fits to rings doesn't quite get the center right, so can shift it manually this way). Both original and new center coordinates are added to the figure's UserData.
%       Press 't' and click on any of the time traces to plot a vertical dotted line. This will appear in all three time traces so you can see if any events happened at the same point in time. Time points are added to the figure's UserData.
%       Press 'm' and click any two points in the time traces to measure the mean value of the coordinate between those time points.
%       Press 'v' and click on any two points in the arc length trace to fit a line across those time points and calculate velocity. Velocities and associated fit results are added to the figure's UserData.
%       Press 's' to save the figure with all associated UserData.
%       Press 'x' to exit the interactive mode.


function verciniSPT_v1(fileFilter, interval, pixSz, varargin)

if nargin == 0 % debug mode
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
    param.fixedRadius = 1; % use fixed radius (otherwise use instantaneous 'radius')

    % filters
    param.time_range = [10 Inf]; % [frames] threshold for how long a track needs to last
    param.pairing_distance_threshold = 2; % [um] threshold for how far a track center can be from a septum center
    param.distance_threshold = 0.05; % [um] threshold for end-to-end distance a track needs to go
    
    param.diffusion_model = 'brownian'; % brownian or directional
    param.msd_fit = 'loglog'; % use loglog fit or linear fit?
    param.msd_fit_range = [0 15]; % range of time-steps to fit for MSD data
    param.R2_threshold = 0.95; % threshold for R^2 (for fit of log(MSD) to log(t))
    
    % no filters
    param.distance_threshold = 0.00; % [um] threshold for end-to-end distance a track needs to go
    param.R2_threshold = -Inf; % threshold for R^2 (for fit of log(MSD) to log(t))
else
    %     Zim_file = uigetfile('*.ome.tif', 'FtsZ image file', [path '\']);
    [im_file, im_path] = uigetfile(fileFilter, 'Pbp2B stack');
    [tracks_file, tracks_path] = uigetfile('*.csv', 'Track file', [im_path '\']);
    [septa_file, septa_path] = uigetfile('*.ome.roi.zip.csv', 'Septa file', [im_path '\']);
    
    param.directory = im_path;
    param.Date = im_file(1:6); % assumes first 6 digits are date
    param.file = im_file(8); % assumes 8th character is file
    
    %     Zim = imreadstack([path '\' Zim_file]);
    im = imreadstack([im_path im_file]);
    tracks = xlsread([tracks_path tracks_file]);
    cells = xlsread([septa_path septa_file]);

    param.interval = interval; % [s]
    param.pixSz = pixSz; % [nm/pix]

    % defaults. Can be overridden by user.
    plot_fov = 0;
    param.fixedRadius = 1; % use fixed radius (otherwise use instantaneous 'radius')
    param.psfFWHM = 250; % [nm]
    param.time_range = [10 Inf]; % [frames] threshold for how long a track needs to last
    param.distance_threshold = 0.00; % [um] threshold for end-to-end distance a track needs to go
    param.pairing_distance_threshold = 2; % [um] threshold for how far a track center can be from a septum center
    param.diffusion_model = 'brownian'; % brownian or directional
    param.msd_fit = 'loglog'; % use loglog fit or linear fit?
    param.msd_fit_range = [0 15]; % range of time-steps to fit for MSD data
    param.R2_threshold = -Inf; % threshold for R^2 (for fit of log(MSD) to log(t))
    
    ii = 1; ringFitArg={};
    while ii<=numel(varargin)
        if strcmp(varargin{ii},'PlotFOV')
            plot_fov = 1;
            ii=ii+1;
        elseif strcmp(varargin{ii},'PsfFWHM')
            param.psfFWHM=varargin{ii+1};
            ii=ii+2;
        elseif strcmp(varargin{ii},'TimeRange')
            param.time_range=varargin{ii+1};
            ii=ii+2;
        elseif strcmp(varargin{ii},'DistThresh')
            param.distance_threshold=varargin{ii+1};
            ii=ii+2;
        elseif strcmp(varargin{ii},'PairDistThresh')
            param.pairing_distance_threshold=varargin{ii+1};
            ii=ii+2;
        elseif strcmp(varargin{ii},'DiffusionModel')
            param.diffusion_model=varargin{ii+1};
            ii=ii+2;
        elseif strcmp(varargin{ii},'MSDfit')
            param.msd_fit=varargin{ii+1};
            ii=ii+2;
        elseif strcmp(varargin{ii},'MSDfitRange')
            param.msd_fit_range=varargin{ii+1};
            ii=ii+2;
        elseif strcmp(varargin{ii},'RsquaredThresh')
            param.R2_threshold=varargin{ii+1};
            ii=ii+2;
        else
            ii=ii+1;
        end
    end

end

pix2um = param.pixSz/1000; % [um/pix]

% sort rows (times sometimes are scrambled)
tracks = sortrows(tracks,[2 7]);

% Start FoV plot
if plot_fov
    figure
    ax_xy = gca;
    hold on
    ylabel('y (\mum)')
    xlabel('x (\mum)')
    set(ax_xy, 'YDir', 'reverse')
    xlim([0 66.56])
    ylim([0 66.56])
    axis equal
end

% make cell array of septal (i.e. cell) coordinates

cells = sortrows(cells, [1 2]); % rearrange so first x coordinate is always lowest, second is highest

Cells = {}; cc1 = 1;
for cc = 1:4:length(cells) % data are by fours, so increment by 4

    Cells{cc1}(:,1) = [cells(cc,3); cells(cc+1,3)];
    Cells{cc1}(:,2) = [cells(cc+1,4); cells(cc+2,4)];

    cc1 = cc1+1;
end

% separate tracks

tracknums = unique(tracks(:,2));

dr=[]; subim=zeros(60,60); subim_2b=zeros(60,60,size(im,3));
for ii = 1:length(Cells)
    
    ud.cell.cellnum = ii;
    
    % SET UP FIGURE FOR THIS CELL
    
    ud.figsavename = [param.directory param.Date '_' param.file '_cell' num2str(ii) '_diffuse'];
    fig_tracks = figure('Position',[850, 130, 770, 800], 'FileName', ud.figsavename);
    ud.axes.ax_im = subplot(4,4,1);
    set(ud.axes.ax_im,"YDir","reverse")
    hold on
    axis equal
    
    ud.axes.ax_polar = subplot(4,3,2,polaraxes);
    hold on
    filename_plot = replace([param.Date '_' param.file], '_', '\_');
    title(ud.axes.ax_polar, ['HaloTag-PBP2B tracking data for FoV ' filename_plot ', cell ' num2str(ii)])
    ud.axes.ax_polar.ThetaDir = 'clockwise';
    ud.axes.ax_polar.ThetaTick = 0:90:360;
    
    ud.axes.ax_msd = subplot(4,3,3);
    hold on
    ylabel('MSD (\mum^2)')
    xlabel('\Deltat (s)')
    set(ud.axes.ax_msd,'YAxisLocation','right')

    ud.axes.ax_rho = subplot(4,3,4:6);
    hold on
    ylabel('\rho (nm)')
    ud.axes.ax_dist = subplot(4,3,7:9);
    hold on
    ylabel('Arc length (nm)')
    
    ud.axes.ax_int = subplot(4,3,10:12);
    hold on
    ylabel('Intensity')
    xlabel('Time (s)')
    
    % FIND COORDINATES FOR CENTER OF RING
    
    % create bounding box (need to cut if near border)
    box_x1 = max([Cells{ii}(1,1), 1]); % [pix]
    box_x2 = min([Cells{ii}(2,1), size(im,1)]); % [pix]
    box_y1 = max([Cells{ii}(1,2), 1]); % [pix]
    box_y2 = min([Cells{ii}(2,2), size(im,2)]); % [pix]
    xdim = box_x2 - box_x1; % [pix]
    ydim = box_y2 - box_y1; % [pix]
    
    ud.cell.boxsize = [xdim ydim];
    
    % get individual ring bounded by the box at index ind_box
    subim(1:ydim+1,1:xdim+1,1:size(im,3)) = im(box_y1:box_y2, box_x1:box_x2, :); % use HT-PBP2B image to find ring. need to flip x and y
%     subim(1:ydim+1,1:xdim+1,1) = Zim(box_y1:box_y2, box_x1:box_x2, 1); % use GFP-FtsZ image to find ring
    
    % fit individual ring to model, get center position
    %ringIm = mean(subim,3);
    ringImMax = max(subim,[],3);
    fitRingArg = {'FixedPositionFit', true};
    try
        [fitParAvg, ~, im_bg_sub] = fitRing(ringImMax, param.pixSz, param.psfFWHM, fitRingArg{:});
    catch
        continue
    end
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
        otsuThresh(:,:,kk) = otsu(imBlur(:,:,ii));
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
    Tracks={}; maxrho=[]; ud.phandles.p_polar=[]; ud.phandles.p_rho=[]; ud.phandles.p_dist=[]; ud.phandles.p_msd=[]; ud.phandles.p_msd_fit=[]; ud.cell.track={};
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
        ud.cell.track{tt}.intensity = Tracks{tt}(:,13) - bg(Tracks{tt}(:,7)+1)';
        
        ud.cell.track{tt}.rho = Tracks{tt}(:,22);
        ud.cell.track{tt}.theta = Tracks{tt}(:,23);
        ud.cell.track{tt}.rho_new = Tracks{tt}(:,22); % duplicate in case you shift it later
        ud.cell.track{tt}.theta_new = Tracks{tt}(:,23); % duplicate in case you shift it later
        
        if param.fixedRadius
            Tracks{tt}(:,24) = Tracks{tt}(:,23) .* radius*param.pixSz; % [nm] angular position, fixed R
        else
            Tracks{tt}(:,24) = Tracks{tt}(:,23) .* Tracks{tt}(:,22); % [nm] angular position, instantaneous R
        end
        Tracks{tt}(:,23) = mod(Tracks{tt}(:,23), 2*pi); % [rad] convert from [-pi, pi] to [0, 2*pi]
        Tracks{tt}(:,23) = Tracks{tt}(:,23) * 180/pi; % [deg]
        
        % patch up boundary conditions
        param.derivative_threshold = 1000 / param.interval; % set threshold based on interval
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
        
        
        % CALCULATE JUMP DISTANCES FOR EACH TIMESTEP (both 1D and 2D)
        
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

            Dt_fit = Dt_cut; % for fitting, use regular Dt (not log)

            [fitvals, res, J] = nlinfit(Dt_fit, MSD_fit, eqn_diff, init_diff);
            ci = nlparci(fitvals, res, 'jacobian', J);
            t = tinv(1-0.05/2, length(MSD_fit)-length(fitvals));
            fiterrs = (ci(:,2)-ci(:,1)) ./ (2*t);
            
            Dt_plot = Dt; % fot plotting the model, use regular Dt (not log)
        end

        sse = sum((log(MSD)-eqn_diff(fitvals,log(Dt))).^2); % residual sum of squares
        sst = sum((log(MSD)-mean(log(MSD))).^2); % total sum of squares
        r2 = 1 - sse/sst; % R squared
        
        if r2 < param.R2_threshold % filter by goodness of fit to MSD
            continue
        end
        
        ud.cell.track{tt}.diffusion_eqn = eqn_diff;
        ud.cell.track{tt}.msd_fits = fitvals; % [D alpha] or [D v], depending on which model you chose
        ud.cell.track{tt}.msd_errs = fiterrs; % [D alpha] or [D v], depending on which model you choose
        
        
        % PLOT EACH TRACK IN POLAR COORDINATES, AS TIME-TRACE, AND ALSO MSD
        
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
    
    set(fig_tracks,'WindowKeyPressFcn',@multikey_analysis_v1)

end

