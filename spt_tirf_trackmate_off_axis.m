% Author: Kevin Whitley
% Date created: 231213

% This function measures the off-axis motion of single-molecule tracks
% along a defined axis. The video file has a track file from TrackMate (in
% csv format) and a septum coordinate file (in roi.zip.csv format). It then
% rotates the track to be along the septal axis and calculates how far it
% deviates from this axis by measuring the distances of each localization
% in the track from this axis.

% The first input into this function is a file. The function assumes there
% will be three associated files for this FoV. Each must have the same
% name, but with different extensions:
%   - The video file (.tif)
%   - The track file outputted from TrackMate (_tracks.csv)
%   - The septa file containing septal coordinates (.roi.zip.csv)
%
% INPUTS:
%   - file: FoV to be analyzed
%   - interval: framerate [s]
%   - pixSz: conversion factor between pixels and nm [nm/pix]
%   - thr_length: threshold for how long a track needs to last to
%       be included [frames]. Default: [0 Inf]
%   - thr_dist_septa: minimum distance a track center can be from a septum
%       to be included [um]. Default: 0.2.
%   - thr_e2e: minimum end-to-end distance a track needs to go to be
%       included [um]. Default: 0.
%   - thr_angle_septa: Minimum angle a track can form with the septal axis
%       to be included [deg]. Default: 30.
%   - plot_xy: plot tracks and septal coordinates in Cartesian coordinates?
%       0: no, 1: yes. Default: 0.
%   - plot_xy_rot: plot tracks and septal coordinates rotated so that each
%       track is plotted along its septal axis. 0: no, 1: yes. Default: 0
% 
% OUTPUT:
%   - OffAxisDat: A structure variable containing information on the FoV.
%       Fields:
%       - datapath: name of the directory the files came from
%       - datafile: name of the file containing tracking data
%       - septafile: name of the file containing septal coordinates
%       - track_dat: all data on tracks in this FoV, pulled from the
%           tracking data file
%       - septa_dat: all data on septa in this FoV, pulled from the septa
%           data file
%       - param: various parameters involved in the analysis (e.g.
%           interval, pixSz, thresholds)
%       - TrackDat: a cell array with information on each track in the FoV:
%           - track_data: data on this track from tracking data file
%           - size: length of track
%           - dist: end-to-end distance traveled by track
%           - min_septa_dist: distance from track to the nearest septum
%           - septum_data: data on the nearest septum to this track
%           - var_cell_ax: variance of the track along the cell long axis
%           - track_septum_angle: angle between a line fitted to the track
%               and the associated septum
%           - wobble_off: mean of distances between each point in the track
%               and the septal axis (mean of rotated tracks)
%           - wobble_fit: mean of distances between each point in the track
%               and the septal axis (fitted line to rotated tracks)


function OffAxisDat = spt_tirf_trackmate_off_axis(file, interval, pixSz, thr_length, thr_dist_septa, thr_e2e, thr_angle_septa, plot_xy, plot_xy_rot)

param.interval = interval; % [s] framerate
param.pix2um = pixSz / 1e3; % [um/pix] pixel to um conversion factor

% FILTERS

if nargin < 4 % default filters
    thr_length = [0 Inf];
    thr_dist_septa = 0.2;
    thr_e2e = 0;
    thr_angle_septa = 30;
end

param.thr_length = thr_length; % [frames] threshold for how long a track needs to last
param.thr_dist_septa = thr_dist_septa; % [um] threshold for how far a track center can be from a septum to be included
param.thr_e2e = thr_e2e; % [um] threshold for end-to-end distance a track needs to go
param.thr_angle_septa = thr_angle_septa; % [deg] threshold for how much a track direction can deviate from the septal axis


if nargin < 8 % default: no plots
    plot_xy = 0;
    plot_xy_rot = 0;
end

track_file = replace(file, '.tif', '_tracks.csv');
septa_file = replace(file, '.tif', '.roi.zip.csv');

if isfile(track_file)
    dat = xlsread(track_file);
else
    error('No track file found. Check file name and extension. (must end with _tracks.csv)')
end
if isfile(septa_file)
    septa = readtable(septa_file);
else
    error('No septa file found. Check file name and extension. (must end with .roi.zip.csv)')
end

%% START PLOTS

if plot_xy
    figure
    p_xy = gca;
    hold on
    ylabel('y (\mum)')
    xlabel('x (\mum)')
    set(p_xy,"YDir","reverse")
end

if plot_xy_rot
    figure
    p_xy_rot = gca;
    hold on
    ylabel('Long axis (\mum)')
    xlabel('Septal axis (\mum)')
end

%% PREPARE COORDINATES OF CELL AXES (LONG AND SHORT)

septa = sortrows(septa, {'Track', 'x'}); % coordinates of septa, output from Export_Coords_to_csv.ijm.
septa.x = septa.x .* param.pix2um; % [um]
septa.y = septa.y .* param.pix2um; % [um]

%% ANALYZE EACH TRACK

% sort rows (times sometimes are scrambled)
dat = sortrows(dat,[2 7]);

tracknums = unique(dat(:,2));

Tracks={}; TrackDat={}; val=[]; JD_all=[]; Septal_all=[]; val_s=[];
for ii = 1:length(tracknums)
    
    Tracks{ii} = dat(dat(:,2)==tracknums(ii),:);
    
    TrackDat{ii}.track_data = dat(dat(:,2)==tracknums(ii),:);
    TrackDat{ii}.size = size(Tracks{ii},1);
    
    % REMOVE TRACKS THAT ARE TOO SHORT OR TOO LONG
    if size(Tracks{ii},1) < param.thr_length(1) || size(Tracks{ii},1) > param.thr_length(2)
        continue
    end
    
    % REMOVE TRACKS THAT DIDN'T GO FAR
    
    Rxy = Tracks{ii}(:,4:5); % [um] (x,y) pairs
    
    TrackDat{ii}.dist = norm(Rxy(end,:)-Rxy(1,:));
    
    if norm(Rxy(end,:)-Rxy(1,:)) < param.thr_e2e
        continue
    end
    
    % FILTER OUT TRACKS THAT AREN'T SEPTAL
    
    % calculate distances of track center to cell long axes
    track_center = [mean(Tracks{ii}(:,4)) mean(Tracks{ii}(:,5))]; % [um]
    
    % calculate distances of track center to septa
    septa_xy1 = [septa.x(1:2:end-1) septa.y(1:2:end-1)]; % [um]
    septa_xy2 = [septa.x(2:2:end) septa.y(2:2:end)]; % [um]
    septa_L = vecnorm(septa_xy1-septa_xy2,2,2); % [um] length of septa
    
    r1_septa = sqrt((track_center(1)-septa_xy1(:,1)).^2 + (track_center(2)-septa_xy1(:,2)).^2); % [um]
    r2_septa = sqrt((track_center(1)-septa_xy2(:,1)).^2 + (track_center(2)-septa_xy2(:,2)).^2); % [um]
    
    th1_septa = acos((r1_septa.^2+septa_L.^2-r2_septa.^2) ./ (2*r1_septa.*septa_L)); % angle between septal axis and r1. law of cosines.
    th2_septa = acos((r2_septa.^2+septa_L.^2-r1_septa.^2) ./ (2*r2_septa.*septa_L)); % angle between septal axis and r1. law of cosines.
    
    th1_septa(th1_septa>pi/2 | th2_septa>pi/2) = nan; % remove any with obtuse angles (otherwise tracks would be found far away from septum, as long as sin(th1)~0)
    th2_septa(th1_septa>pi/2 | th2_septa>pi/2) = nan; % remove any with obtuse angles (otherwise tracks would be found far away from septum, as long as sin(th1)~0)
    
    dist1_septa = r1_septa .* sin(th1_septa); % shortest distance of track center to all septa
    dist2_septa = r2_septa .* sin(th2_septa); % shortest distance of track center to all septa (same as dist1_septa - but good to check)
    
    TrackDat{ii}.min_septa_dist = min(dist1_septa);
    
    % remove tracks not near septa
    [val_s(ii), ind_s] = min(dist1_septa); % minimum distance between track and septa
    if val_s(ii) > param.thr_dist_septa || isnan(val_s(ii))
        continue
    end
    
    % CALCULATE VARIANCE ACROSS CELL LONG AXIS
    
    % rotate track coordinates to be along septal axis
    s_vec = septa_xy2(ind_s,:) - septa_xy1(ind_s,:); % vector of septum
    s_th = atan(s_vec(2)/s_vec(1)); % [rad] angle of septum relative to (x,y) system
    s_th = s_th*180/pi; % [deg]
    R = rotz(-s_th); % rotation matrix
    
    track_xyz = [Rxy'; zeros(1,length(Rxy))]; % [um]
    track_rot = R * track_xyz; % [um] track coordinates rotated to be along septal axis
    
    % get variance along cell long axis
    
    var_cell_ax = var(track_rot(2,:)); % [um] variance along cell long axis
    
    TrackDat{ii}.septum_data = [septa_xy1(ind_s,:); septa_xy2(ind_s,:)];
    TrackDat{ii}.var_cell_ax = var(track_rot(2,:));
    
    % fit line to track
    lineqn = @(a,x) a(1) + a(2)*x;
    [m, res] = nlinfit(track_rot(1,:), track_rot(2,:), lineqn, [0 0]);
    
    
    % remove tracks with angle deviating too much from septal axis (can
    % happen if an immobile molecule is tracked along the cell long axis as
    % the cell elongates)
    track_fit_vec = [1 m(2)]; % vector of line fitted to track
    track_sep_dot = dot([1 0], track_fit_vec); % dot product with septal axis (now x axis)
    track_sep_th = acos(track_sep_dot/norm(track_fit_vec)); % [rad] angle between them
    track_sep_th = track_sep_th * 180/pi; % [deg]
    
    TrackDat{ii}.track_septum_angle = track_sep_th;
    
    if track_sep_th > param.thr_angle_septa
        continue
    end
    
    
    % measure deviation of track from its fitted line ('wobble')
    wobbleoff = mean(abs(track_rot(2,:)-mean(track_rot(2,:)))); % [um]
    wobble = mean(abs(res)); % [um]
    
    TrackDat{ii}.wobble_off = wobbleoff; % residuals from offset
    TrackDat{ii}.wobble_fit = mean(abs(res)); % residuals from line fit
    
    % PLOT EVERYTHING
    
    if plot_xy
        % plot track and associated septum
        plot(p_xy, Tracks{ii}(:,4), Tracks{ii}(:,5)) % track
        plot(p_xy, [septa_xy1(ind_s,1) septa_xy2(ind_s,1)], [septa_xy1(ind_s,2) septa_xy2(ind_s,2)]) % septum
    end
    if plot_xy_rot
        this_sep = [[septa_xy1(ind_s,:)'; 1] [septa_xy2(ind_s,:)'; 1]];
        sep_rot = R*this_sep; % rotated septum
        
        plot(p_xy_rot, track_rot(1,:), track_rot(2,:))
        plot(p_xy_rot, [sep_rot(1,2) sep_rot(1,1)], [sep_rot(2,2) sep_rot(2,1)])
        
        tv = min([sep_rot(1,2) sep_rot(1,1)]):0.05:max([sep_rot(1,2) sep_rot(1,1)]);
        plot(p_xy_rot, tv, lineqn(m,tv), 'r')
        
        plot(p_xy_rot, tv, mean(track_rot(2,:)).*ones(1,length(tv)), 'c')
    end
end

OffAxisDat.datapath = path;
OffAxisDat.datafile = track_file;
OffAxisDat.septafile = septa_file;
OffAxisDat.track_dat = dat;
OffAxisDat.septa_dat = septa;
OffAxisDat.param = param;
OffAxisDat.TrackDat = TrackDat;


