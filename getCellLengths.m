% Author: Kevin Whitley
%
% This script measures the lengths of cells from intensity profiles running
% lengthwise along them. To do this, it simply locates the two highest
% peaks along a given line profile and measures the distance between them.
% 
% Data exists in csv files with the format:
% - Odd columns: Distance [um]
% - Even columns: Intensities
% 
% Each odd column and even column exist as pairs (intensities are
% associated with the distances).
% 
% INPUT:
%  - path: Directory where csv files are located.
%  - plotdat: Plot data? 0 = no, 1 = yes.
% 
% OUTPUT:
%  - all_lengths: Array of cell lengths [um].

function all_lengths = getCellLengths(path, plotdat)

if nargin == 0
    path = uigetdir('\\campus\rdw\FMS CBCB\nkw81\');
    plotdat = 1;
end

if nargin == 1
    plotdat = 0;
end

if plotdat
    figure
    xlabel('Distance (\mum)')
    ylabel('Intensity (arb. units)')
    hold on
    box on
end

csvfiles = dir([path, '\*.csv']);
ncsv = length(csvfiles);

all_lengths = [];
for cc = 1:ncsv
    dat = readmatrix([path '\' csvfiles(cc).name]);
    
    leng = [];
    for ii = 1:2:size(dat,2)
        linepro = dat(:,ii+1);
        linepro = (linepro-min(linepro)) / max(linepro);
        [pks,locs,w,p] = findpeaks(linepro);
        pkmat = [pks locs dat(locs,ii) w p];
        pkmat = sortrows(pkmat, 5, 'descend');
        leng(ii) = abs(pkmat(2,3) - pkmat(1,3)); % cell length
        
        if plotdat
            cla
            plot(dat(:,ii),linepro)
            hold on
            plot(pkmat(2,3), pkmat(2,1), '+m')
            plot(pkmat(1,3), pkmat(1,1), 'xm')
        end
    end
    leng(leng==0)=[];
    all_lengths = [all_lengths leng];
end