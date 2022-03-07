
ud = get(gcf,'UserData');

for tt = 1:length(ud.cell.track)
    
    Date = ud.param.Date;
    filename = ud.param.file;
    cellnum = ud.cell.cellnum;
    trackind = ud.cell.track{tt}.trackind;
    tracknum = ud.cell.track{tt}.tracknum;
    
    if isfield(ud.cell.track{tt},'t_pts') && length(ud.cell.track{tt}.t_pts)>1
        t_run = [t_run ud.cell.track{tt}.t_pts(2:end)-ud.cell.track{tt}.t_pts(1:end-1)];
    end
    
    if isfield(ud.cell.track{tt},'speeds')
        spds_all = [spds_all; ud.cell.track{tt}.speeds];
%         runs = size(ud.cell.track{tt}.speeds,1);
        
%         for rr = 1:runs
%             res_all = {Date filename cellnum trackind tracknum t_run rr ud.cell.track{tt}.speeds(rr,:)};
%             
%             Results = {Results; res_all};
%         end
    end
    
end

% var_names = {'Date','FileName','CellNum','TrackIndex','TrackNum','t_run','RunNum','d_0','d_end','t_0','t_end','speed'};