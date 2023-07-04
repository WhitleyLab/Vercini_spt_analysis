% Author: Kevin Whitley
% Date created: 220614

% This script takes the data out of a structure variable containing
% single-particle tracking data, performs analysis, and plots results.

dat = phmm30;

% unpack some relevant data in structure variable and repack it into a
% table so we can work with it more easily.

Results=[];
for ii = 1:size(dat,1)
    
    Date = str2double(dat(ii).param.Date);
    File = str2double(dat(ii).param.file);
    AnalysisDate = str2double(dat(ii).param.analysisDate);
    pixSz = dat(ii).param.pixSz; % [nm/pix]
    Interval = dat(ii).param.interval; % [s]
    
    CellNum = dat(ii).cell.cellnum;
    Radius = dat(ii).cell.cell_radius_fit; % [pix]
    Ntracks = dat(ii).cell.Ntracks;
    
    for jj = 1:size(dat(ii).cell.track,2)
        
        if ~dat(ii).cell.track{jj}.isplotted % didn't do any analysis on this track - ignore
            continue
        end
        
        TrackInd = dat(ii).cell.track{jj}.trackind;
        TrackLength = length(dat(ii).cell.track{jj}.time); % [frames]
        D = dat(ii).cell.track{jj}.msd_fits(1); % [um^2/s]
        vMSD = dat(ii).cell.track{jj}.msd_fits(2); % [um/s]
        
        if isfield(dat(ii).cell.track{jj},'segment')
            
            Nsegments = size(dat(ii).cell.track{jj}.segment,2);
            
            Tdist=0; Cdist=0; Tlife=0;
            for kk = 1:Nsegments
                
                SegNum = kk;
                Pos1 = dat(ii).cell.track{jj}.segment{kk}.pos1; % [nm]
                Pos2 = dat(ii).cell.track{jj}.segment{kk}.pos2; % [nm]
                Distance = Pos2 - Pos1; % [nm]
                Tdist = Tdist + Distance; % [nm]
                Cdist = Cdist + abs(Distance); % [nm]
                
                t1 = dat(ii).cell.track{jj}.segment{kk}.t1; % [s]
                t2 = dat(ii).cell.track{jj}.segment{kk}.t2; % [s]
                Lifetime = t2 - t1; % [s]
                Tlife = Tlife + Lifetime; % [s]
                
                Velocity = dat(ii).cell.track{jj}.segment{kk}.velocity; % [nm/s]
                
                Chi2 = dat(ii).cell.track{jj}.segment{kk}.chi2; % [nm]
                
                if kk == Nsegments
                    TotalDistance = Tdist;
                    CumulativeDistance = Cdist;
                    TotalLifetime = Tlife;
                else
                    TotalDistance = nan;
                    CumulativeDistance = nan;
                    TotalLifetime = nan;
                end
                
                res = table(Date,File,AnalysisDate,pixSz,Interval,CellNum,Radius,Ntracks,...
                    TrackInd,TrackLength,D,vMSD,Nsegments,SegNum,Pos1,Pos2,Distance,t1,t2,...
                    Lifetime,Velocity,Chi2,TotalDistance,CumulativeDistance,TotalLifetime);
                
                Results = [Results; res];
            end
           
        else
            Nsegments = 0;
            SegNum = 0;
            Pos1 = nan;
            Pos2 = nan;
            t1 = nan;
            t2 = nan;
            Velocity = nan;
            Chi2 = nan;
            TotalDistance = nan;
            CumulativeDistance = nan;
            TotalLifetime = nan;
            
            res = table(Date,File,AnalysisDate,pixSz,Interval,CellNum,Radius,Ntracks,...
                    TrackInd,TrackLength,D,vMSD,Nsegments,SegNum,Pos1,Pos2,Distance,t1,t2,...
                    Lifetime,Velocity,Chi2,TotalDistance,CumulativeDistance,TotalLifetime);
                
            Results = [Results; res];
        end
    end

end


