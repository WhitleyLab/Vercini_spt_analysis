% Author: Kevin Whitley
% Date created: 220613

% This function takes a table outputted from
% Batch_measure_kymotraces_polyline.ijm and returns a table with only track
% segments that are 'true'. That is, it removes track segments where:
%       - The intensity dropped to zero (which has a high chance of
%       resulting from photobleaching rather than dissociation)
%       - The segment reached the end of the acquisition interval
%       - The segment began with the acquisition interval

function [Tout, nsegs] = getKymotraceMiddleSegments(Tin, interval)

% first, remove all global measurements (only want track segments, not full
% tracks)
T = Tin(Tin.Segment_>0,:);

% next, go through each trace (each ROI) individually
[~, inds] = unique(T.Image_ROI_Name);
inds = sort(inds);

Tout = []; nsegs = [];
for ii = 1:length(inds)-1
    
    % for each trace (each ROI), go through segments
    segs = T(inds(ii):inds(ii+1)-1,:);
    
    runsegs = segs(abs(segs.Speed_nm_s_)>5 & abs(segs.Speed_nm_s_)<50,:);
    nsegs(ii) = height(runsegs);
    
    % exclude segments that started with acquisition interval (within one
    % frame)
    segs(segs.y1_s_<interval,:) = [];

    % exclude all end segments (may have ended in photobleaching). this
    % also removes any that ended with the acquisition interval.
    if ~isempty(segs)
        segs(segs.Segment_(:)==segs.Segment_(end),:) = [];
    end
    
    Tout = [Tout; segs];
end

