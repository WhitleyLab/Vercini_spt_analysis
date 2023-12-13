% Author: Kevin Whitley
% Date created: 231115
% Updated: 231213

% This script filters horizontal cell SPT data to get only directional
% tracks at septa, and then makes arrays of the variance along the cell
% axis and the linear fit residuals ('wobble').
% 
% The structure variable AllDat used as an input in this script is the
% output from spt_tirf_trackmate_off_axis or
% spt_tirf_trackmate_off_axis_batch2.
%
% Note that filters exist both here and in the code preceding this
% (spt_tirf_trackmate_off_axis_batch2). You can filter the data in the
% previous code, or instead do it here.
% 
% INPUTS:
%   - AllDat: Structure variable outputted from spt_tirf_trackmate_off_axis
%       or spt_tirf_trackmate_off_axis_batch2
%   - thr_length: threshold for how long a track needs to last to
%       be included [frames]. Default: [0 Inf]
%   - thr_dist_septa: minimum distance a track center can be from a septum
%       to be included [um]. Default: 0.2.
%   - thr_e2e: minimum end-to-end distance a track needs to go to be
%       included [um]. Default: 0.
%   - thr_angle_septa: Minimum angle a track can form with the septal axis
%       to be included [deg]. Default: 30.
%   - plot_hist: plot standard deviations and wobbles from analysis here.
%       0: no, 1: yes. Default: 0.
%
% OUTPUTS:
%   - var_cell_ax: vector of variances along the cell long axis [um^2]
%   - wobble_fit: vector of distances from track localizations from line
%       fitted to track [um]
%   - wobble_off: vector of distances from track localizations from septal
%       axis [um]


function [var_cell_ax, wobble_fit, wobble_off] = plot_spt_tirf_trackmate_off_axis2(AllDat, thr_length, thr_dist_septa, thr_e2e, thr_angle_septa, plot_hist)

% FILTERS

if nargin < 2 % default filters
    thr_length = [10 Inf]; % [frames] threshold for how long a track needs to last
    thr_dist_septa = 0.2; % [um] threshold for how far a track center can be from a septum to be included
    thr_e2e = 0.00; % [um] threshold for end-to-end distance a track needs to go
    thr_angle_septa = 30; % [deg] threshold for how much a track direction can deviate from the septal axis
end

if nargin < 6
    plot_hist = 0;
end

param.thr_length = thr_length;
param.thr_dist_septa = thr_dist_septa;
param.thr_e2e = thr_e2e;
param.thr_angle_septa = thr_angle_septa;

%% FILTER DATA AND GET VARIANCE AND WOBBLES

var_cell_ax=[]; wobble_fit=[]; wobble_off=[];
for ff = 1:length(AllDat)
    
    TrackFoV = AllDat(ff); % tracks for one FoV
    
    for ii = 1:length(TrackFoV.TrackDat)
        
        track = TrackFoV.TrackDat{ii}; % single track
        
        % REMOVE TRACKS THAT ARE TOO SHORT OR TOO LONG
        if track.size < param.thr_length(1) || track.size > param.thr_length(2)
            continue
        end
        
        % REMOVE TRACKS THAT DIDN'T GO FAR
        if track.dist < param.thr_e2e
            continue
        end
        
        % FILTER OUT TRACKS THAT AREN'T SEPTAL
        if track.min_septa_dist > param.thr_dist_septa || isnan(track.min_septa_dist)
            continue
        end
        
        % FILTER OUT TRACKS THAT DIDN'T MOVE ALONG SEPTAL AXIS
        if track.track_septum_angle > param.thr_angle_septa
            continue
        end
        
        var_cell_ax = [var_cell_ax; track.var_cell_ax]; % [um^2]
        wobble_fit = [wobble_fit; track.wobble_fit]; % [um]
        wobble_off = [wobble_off; track.wobble_off]; % [um]
        
    end
    
end

%% PLOT HISTOGRAMS OF SD AND WOBBLES

if plot_hist
    
    figure
    histogram(sqrt(var_cell_ax.*1e6), 0:5:100)
    xlabel('Standard deviation along cell axis (nm)')
    ylabel('Counts')
    
    figure
    histogram(wobble_fit.*1e3)
    xlabel('Mean of distances from line fitted to track (nm)')
    ylabel('Counts')
    
    figure
    histogram(wobble_off.*1e3)
    xlabel('Mean of distances from septal axis (nm)')
    ylabel('Counts')
    
end