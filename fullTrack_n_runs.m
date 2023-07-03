% Author: Kevin Whitley
% Date created: 221220
%
% This function takes a table outputted from
% Batch_measure_kymotraces_polyline.ijm and gets the number of run segments
% per full track. It outputs these in a new table.

function T2 = fullTrack_n_runs(T, spdcut)

if nargin < 2
    spdcut = [0 Inf];
end

if ismember('n_segs', T.Properties.VariableNames)
    T = removevars(T, 'n_segs');
end
if ismember('n_runs', T.Properties.VariableNames)
    T = removevars(T, 'n_runs'); % clear out the old column
end

% first, grab a track with all its segments (and global measurements)
[~, inds] = unique(T.("Image_ROI_Name"));
inds = sort(inds);

n_segs = nan(size(T,1),1); n_runs=nan(size(T,1),1);
for ii = 1:length(inds)
    trackname = T{inds(ii),'Image_ROI_Name'}{1};
    track = T(strcmp(T.("Image_ROI_Name"),trackname),:);
    
    track_segments = track(track.Segment_>0,:); % segment measurements only, no global
    
    if any(abs(track_segments.Speed_nm_s_)>spdcut(2)) % if there are any diffusive segments, don't use this track
        n_segs(inds(ii)) = nan;
        n_runs(inds(ii)) = nan;
    else
        n_segs(inds(ii)) = height(track_segments);
        n_runs(inds(ii)) = sum(abs(track_segments.Speed_nm_s_)>spdcut(1) & abs(track_segments.Speed_nm_s_)<spdcut(2)); % total number of runs in this track
    end
end

T2 = addvars(T, n_segs);
T2 = addvars(T2, n_runs);