% Author: Kevin Whitley
% Date created: 231109

% This function measures the off-axis motion of single-molecule tracks
% along a defined axis. Each video file has a track file from TrackMate (in
% csv format) and a septum coordinate file (in roi.zip.csv format). It then
% rotates the track to be along the septal axis and calculates how far it
% deviates from this axis by measuring the distances of each localization
% in the track from this axis.

% The first input into this function is a directory. The function assumes
% there will be three associated files for each FoV you'll analyze. Each
% must have the same name, but with different extensions:
%   - The video file (.tif)
%   - The track file outputted from TrackMate (_tracks.csv)
%   - The septa file containing septal coordinates (.roi.zip.csv)
%
% INPUTS:
%   - path: path to a directory containing the files to analyze in batch
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
%   - AllDat: A structure variable containing information on each FoV.
%       Each FoV has fields:
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

function AllDat = spt_tirf_trackmate_off_axis_batch2(path, interval, pixSz, thr_length, thr_dist_septa, thr_e2e, thr_angle_septa, plot_xy, plot_xy_rot)

% FILTERS

% if nargin < 4 % default filters
%     thr_length = [0 Inf];
%     thr_dist_septa = 0.2;
%     thr_e2e = 0;
%     thr_angle_septa = 30;
% end
% 
% if nargin < 8 % default: no plots
%     plot_xy = 0;
%     plot_xy_rot = 0;
% end

vid_files = dir([path '\*.tif']);

AllDat=[];
for fnum = 1:length(vid_files)
    
    track_file = replace(vid_files(fnum).name, '.tif', '_tracks.csv');
    septa_file = replace(vid_files(fnum).name, '.tif', '.roi.zip.csv');
    
    if isfile([path '\' track_file])
        dat = xlsread([path '\' track_file]);
    else
        continue
    end
    if isfile([path '\' septa_file])
        septa = readtable([path '\' septa_file]);
    else
        continue
    end
    
    OffAxDat = spt_tirf_trackmate_off_axis([vid_files(fnum).folder '\' vid_files(fnum).name], interval, pixSz, thr_length, thr_dist_septa, thr_e2e, thr_angle_septa, plot_xy, plot_xy_rot);
    
    AllDat(fnum).datapath = OffAxDat.datapath;
    AllDat(fnum).datafile = OffAxDat.datafile;
    AllDat(fnum).septafile = OffAxDat.septafile;
    AllDat(fnum).track_dat = OffAxDat.track_dat;
    AllDat(fnum).septa_dat = OffAxDat.septa_dat;
    AllDat(fnum).param = OffAxDat.param;
    AllDat(fnum).TrackDat = OffAxDat.TrackDat;

end
