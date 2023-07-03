% Date: 230626

% This script will grab all jump-distance data from files in 220503 and
% 220505 folders (30 ms acquisition vercini data) and put them into a
% matrix of 300 columns.

% JD_mat = NaN(1,300); % initialize matrix

path = uigetdir('\\campus\rdw\FMS CBCB\nsh167\Shared\data\Whitley-Kevin\');
figfiles = dir([path '\*cell*.fig']);

for ii = 1:length(figfiles)
    
    open([figfiles(ii).folder '\' figfiles(ii).name])
    
    ud = get(gcf,'UserData');
    
    for jj = 1:length(ud.cell.track)
        if isfield(ud.cell.track{jj},'jump_distances')
            udmat = [ud.cell.track{1,jj}.jump_distances NaN(size(ud.cell.track{1,jj}.jump_distances,1), 300-size(ud.cell.track{1,jj}.jump_distances,2))];
            
            JD_mat = [JD_mat; udmat];
        elseif isfield(ud.cell.track{jj},'jump_distances_2D') % 220505 data
            udmat = [ud.cell.track{1,jj}.jump_distances_2D NaN(size(ud.cell.track{1,jj}.jump_distances_2D,1), 300-size(ud.cell.track{1,jj}.jump_distances_2D,2))];
            
            JD_mat = [JD_mat; udmat];
        end
    end

    close
    
end