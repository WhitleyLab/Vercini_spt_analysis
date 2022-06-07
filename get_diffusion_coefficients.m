
ud = get(gcf,'UserData');

for tt = 1:length(ud.cell.track)
    
    Date = ud.param.Date;
    filename = ud.param.file;
    cellnum = ud.cell.cellnum;
    trackind = ud.cell.track{tt}.trackind;
    tracknum = ud.cell.track{tt}.tracknum;
    
    if ud.cell.track{tt}.isplotted
        D = [D exp(ud.cell.track{tt}.msd_fits(1))];
        alpha = [alpha ud.cell.track{tt}.msd_fits(2)];
    end
    
end

% var_names = {'Date','FileName','CellNum','TrackIndex','TrackNum','t_run','RunNum','d_0','d_end','t_0','t_end','speed'};