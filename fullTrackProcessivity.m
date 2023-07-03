% Author: Kevin Whitley
% Date created: 220606
%
% This function takes a table outputted from
% Batch_measure_kymotraces_polyline.ijm and gets the full track
% processivities. It outputs these in a new table.

function T2 = fullTrackProcessivity(T, spdcut)

if nargin < 2
    spdcut = [0 Inf];
end

if ismember('Total_Length_nm_', T.Properties.VariableNames)
    T = removevars(T, 'Total_Length_nm_'); % clear out the old column
end

% first, grab a track with all its segments (and global measurements)
[~, inds] = unique(T.("Image_ROI_Name"));
inds = sort(inds);

Total_Length_nm_=nan(size(T,1),1);
for ii = 1:length(inds)
    trackname = T{inds(ii),'Image_ROI_Name'}{1};
    track = T(strcmp(T.("Image_ROI_Name"),trackname),:);
    
    track_segments = track(track.Segment_>0,:); % segment measurements only, no global
    
    if ~any(abs(track_segments.Speed_nm_s_)<spdcut(1)) && ~any(abs(track_segments.Speed_nm_s_)>spdcut(2))
        Total_Length_nm_(inds(ii)) = sum(abs(track_segments.Length_nm_));
    end
end

T2 = addvars(T, Total_Length_nm_);