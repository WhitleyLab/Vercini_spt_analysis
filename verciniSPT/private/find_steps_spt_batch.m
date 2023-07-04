function find_steps_spt_batch
% find state changes in single particle tracking trajectories using
% AutoStepFinder on the derivative of the position trace (assuming discrete
% changes in velocity). this version does all trajectories in batch.

ud = get(gcf,'UserData');

for ii = 1:length(ud.phandles.p_dist)
    
    if ~ud.cell.track{ii}.isplotted
        continue
    end

    time = get(ud.phandles.p_dist(ii),'xdata');
    dist = get(ud.phandles.p_dist(ii),'ydata');
    
    % get velocity vs. time (first derivative of distance)
    % dsm = smooth(dist, 2)';
    dsfact = 1; % for downsampling
    % dsm = downsample(dist, dsfact)';
    % tsm = downsample(time, dsfact)';
    
    wdw = 1;
    vel = (dist(2:end)-dist(1:end-1)) ./ (time(2:end)-time(1:end-1));
    
    % use 5-point stencil to approximate derivative (less noisy)
    % window = 5;
    % vel = (-dist(5:end) + 8*dist(4:end-1) - 8*dist(2:end-3) + dist(1:end-4)) / 12;
    
    % use savitzky-golay filter
%     order = 2;
%     wdw = 3;
%     [~, dfilt] = sgolay(order,wdw);
%     vel = conv(dist, factorial(1)/(-1)^1 * dfilt(:,2), 'same');
%     vel = vel(wdw-1:end-(wdw-2));
    % tsm = time;
    
    
    % vel = (dsm(2:end)-dsm(1:end-1)) ./ (tsm(2:end)-tsm(1:end-1)); % [nm/s]
    % vel = (dist(2:end)-dist(1:end-1)) ./ (time(2:end)-time(1:end-1)); % [nm/s]
    
    
    
    % write velocity coordinates to temporary txt file for AutoStepFinder to use
    fileid = fopen([ud.figsavename '.txt'], 'w');
    fprintf(fileid, '%d\r\n',vel);
    fclose(fileid);
    
    % run AutoStepFinder (without GUI)
    init.datapath = ud.param.directory;
    init.txtfile = [ud.figsavename '.txt'];
    init.codefolder = 'C:\Users\nkw81\Documents\GitHub\Vercini_spt_analysis\';
    init.SMaxTreshold = 0.15;
    init.overshoot       = 1;
    init.fitrange = 10000;
    init.stepnumber = 10000;
    init.resolution = ud.param.interval * dsfact;
    init.nextfile = 1;
    init.meanbase = 0;
    init.max_range = 100;
    init.userplt         = 0;           %Turn user plot function on/ off
    init.scurve_eval     = 1;          %Turn S-curve evaluation on/ off
    init.fitmean         = 1;             %Use mean for fitting
    init.fitmedian       = 0;           %Use median for fitting
    init.treshonoff      = 0;         %Turn base line treshholding on/ off
    init.txtoutput       = 0;           %Output .txt files
    init.matoutput       = 1;           %Output .mat files
    init.setsteps        = 10; %Resolution of measurement
    init.parametersout   = 0;       %Save parameters output file.
    init.fitsoutput      = 1;          %Save fits output file.
    init.propoutput      = 1;          %Save properties output file.
    init.scurvesoutput   = 1;       %Save S-curves output file.
    init.manualoff       = 1;           %Manual mode off
    init.manualon        = 0;            %Manual mode on
    init.estimatenoise   = 0;         %Noise estimation on
    if init.treshonoff   == 1
        init.basetresh     = init.meanbase;                         %Treshhold the mean of your base line
    else
        init.basetresh     = -100000;
    end
    
    init.booton   = 0;         %Noise estimation on
    if init.booton   == 1
        init.bootstraprepeats=1000;                           %Add bootstrap erorrs per step
    else
        init.bootstraprepeats=0;
    end
    init.singlerun       = 1;             %Single or batch run
    if init.singlerun    == 1
        init.hand_load     =  1;                                       %Single Run
        init.rerun         =  0;
        if init.rerun      == 1
            init.hand_load =  0;
        end
    else
        init.hand_load     =  2;                                       %Batch Run
        init.datapath      = uigetdir(init.datapath);               %Get directory for batch analysis
    end
    
    mainplot = gcf;
    
    stepplot = figure;
    hand.plot_fit = subplot(3, 2, 1:2);
    hand.plot_Scurve = subplot(3, 2, [3 5]);
    hand.plot_results = subplot(3, 2, [4 6]);
    
    try
        autostepfind_nogui(init, hand);
    catch
        % delete temporary txt file
        delete(init.txtfile)
        
        % delete step figure
        close(stepplot)
        
        continue
    end
    
    % delete temporary txt file
    delete(init.txtfile)
    
    % delete step figure
    close(stepplot)
    
    % grab output files and add to userdata
    if exist([ud.param.directory 'StepFit_Result\'], 'dir')
        fit_result = open([ud.param.directory 'StepFit_Result\' ud.param.Date '_' ud.param.file '_cell' num2str(ud.cell.cellnum) '_fits.mat']);
        step_out = open([ud.param.directory 'StepFit_Result\' ud.param.Date '_' ud.param.file '_cell' num2str(ud.cell.cellnum) '_properties.mat']);
        s_curve = open([ud.param.directory 'StepFit_Result\' ud.param.Date '_' ud.param.file '_cell' num2str(ud.cell.cellnum) '_s_curve.mat']);
        
        % delete output folder
        rmdir([ud.param.directory 'StepFit_Result\'], 's');
    else
        fit_result.Time = time;
        fit_result.Data = vel;
        
        X = [time' ones(length(time),1)];
        b = X \ dist'; % linear regression
        
        step_out.IndexStep = [];
        step_out.TimeStep = [];
        step_out.LevelBefore = [];
        step_out.LevelAfter = b(1);
        step_out.StepSize = [];
        step_out.DwellTimeStepBefore = [];
        step_out.DwellTimeStepAfter = [];
        step_out.StepError = [];
        
        s_curve = [];
    end
    
    ud.cell.track{ii}.stepfit.window = wdw;
    ud.cell.track{ii}.stepfit.fit_result = fit_result;
    ud.cell.track{ii}.stepfit.step_output = step_out;
    ud.cell.track{ii}.stepfit.s_curve = s_curve;
    
    % plot result on top of track
    plot_line_segment_ASF(time, dist, ud.param.interval, step_out, ud.axes.ax_dist, wdw)
    
end


set(mainplot,'userData',ud)