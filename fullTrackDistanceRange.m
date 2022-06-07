% Author: Kevin Whitley
% Date created: 220607
%
% This function takes a table outputted from
% Batch_measure_kymotraces_polyline.ijm and gets the full track distance
% range (minimum position to maximum position). It outputs these in a new
% table.

function T2 = fullTrackDistanceRange(T, spdcut)

if nargin < 2
    spdcut = [0 Inf];
end

% first, grab a track with all its segments (and global measurements)
[~, inds] = unique(T.("Image_ROI_Name"));
inds = sort(inds);

Max_Track_Dist_nm_=nan(size(T,1),1);
for ii = 1:length(inds)
    trackname = T{inds(ii),'Image_ROI_Name'}{1};
    track = T(strcmp(T.("Image_ROI_Name"),trackname),:);
    
    track_segments = track(track.Segment_>0,:); % segment measurements only, no global
    
    if ~any(abs(track_segments.Speed_nm_s_)<spdcut(1)) && ~any(abs(track_segments.Speed_nm_s_)>spdcut(2))
        min_dist = min(track_segments.Length_nm_);
        min_pos = min([0 min_dist]); % position furthest to 'left' (most negative)
        max_dist = max(track_segments.Length_nm_);
        max_pos = max([0 max_dist]); % position furthest to 'right' (most positive)
        
        Max_Track_Dist_nm_(inds(ii)) = max_pos - min_pos;
    end
end

T2 = addvars(T, Max_Track_Dist_nm_);