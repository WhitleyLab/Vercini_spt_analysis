function [ ringStack_noBg, ringKymograph, circleData, kymoInfo, rawKymograph, fitPar] = doBgSubAndKymo(ringStack,pixSzNm,lineWidthNm, psfFWHM,varargin)
% function [ ringStack_noBg, ringKymograph, circleData, kymoInfo] = doBgSubAndKymo(ringStack,pixSzNm,lineWidthNm, psfFWHM,varargin)
%extract circular kymograph, integrating over widthNM annulus thickness
%use fitting to the ring to find the diameter, should make it more robust eg on small rings
%   Performs background subtraction and kymograph fitting for vertically immobilized cells
%   Kymographs are background subtracted ie zero should equal a genuine gap in the ring
% INPUTS:
%   ringStack: [y,x,frame] ring dynamics movie 
%   pixSz: Camera pixel size in nanometres
%   lineWidthNm: perpendicular distance over which to integrate line profile signal to improve SNR. 
%   psfFWHM: Fitted PSF FWHM, nm. Determines PSF size used to blur the fitted ring. 
% 
% OPTIONAL INPUTS:
%   'PsfWidthRangeNm',psfWidthExtraNm: Wiggle room allowed on fitted PSF FHWM. Ie fitted PSF width can be within range psfFWHM +/- psfWidthExtraNm. DEFAULT: 50
%   'CytoplasmBG-FWHM', cytoBgFWHM_nm: Initial guess for the FWHM of the large gaussian fitted to account for the defocussed cytoplasmic background. DEFAULT:1300
%   'CytoplasmBG-FWHM-min', cytoBgFWHMmin_nm: Minimum for the FWHM of the large gaussian fitted to account for the defocussed cytoplasmic background. DEFAULT:1000
%   'CytoplasmBG-FWHM-max', cytoBgFWHMmax_nm: Maximum for the FWHM of the large gaussian fitted to account for the defocussed cytoplasmic background. DEFAULT:1000
% NOTE: A second defocussed Gaussian is also fitted, with min width cytoBgFWHMmin_nm, and max width=Inf because a single gaussian does not fit well the cytoplasmic BG distribution.
%   'RingRadius-max', radMax_nm: Maximum fitted ring radius. Default should hold well for WT or even most mutant Bsubtilis but change if in a different organism. If you set it too large the fitting becomes unstable for small rings. DEFAULT: 600
%   'ZeroPadKymograph', doZeroPadKymo: Add a zero row as the last row of the kymograph so that ImageJ plotting defaults to the correct contrast. DEFAULT:false 
%   'FixedRadiusFit', fitParAvg: Fix the ring radius and shape parameters to the average ring parameters. Note the ring centroid can still shift to allow for small drifts. Useful for cells that dont constrict within timeframe of imaging. If the cells constrict you need to turn this off. FITPARAVG is the result of a prior fit to an averaged ring, used to fix the positions. DEFAULT: true
%   'FitMaxIP', true/false: Use maximum intensity projection instead of average for the fixed radius/ position fitting
%   'HoughCircleGuess',true/false:Uses hough circle finding estimator for the initial guess. In general seems more robust than the simple binary estimator
%
% NOTE: If the background subtration fails for some frames - slow fitting, bright bands in the kymographs - this is usually because 'FixedRadiusFit' is set to true, but the radius is changing  - try changing 'FixedRadiusFit' to false. If the radius is changign the FixedRadiusFit gives bad results as the average radius is not a good match for all frames

nargin = numel(varargin);
fitRingArg={};
doZeroPadKymo = false;
doFixedRadiusFit = true;
doMaxFit=false;
ii = 1;
while ii<=numel(varargin)
    if strcmp(varargin{ii},'ZeroPadKymograph')
        doZeroPadKymo=varargin{ii+1};
        ii=ii+2;
    elseif strcmp(varargin{ii},'FixedRadiusFit')
        doFixedRadiusFit=varargin{ii+1};
        ii=ii+2;
    elseif strcmp(varargin{ii},'FitMaxIP')
        doMaxFit=varargin{ii+1};
        ii=ii+2;
    else
        fitRingArg={fitRingArg{:},varargin{ii}};
        ii=ii+1;
    end
end

ringStack = double(ringStack);

%optionally find the radius etc from the average of whole movie
if doFixedRadiusFit
    if doMaxFit
        ringIm = max(ringStack,[],3);
    else
        ringIm = mean(ringStack,3);
    end
    fitParAvg = fitRing(ringIm, pixSzNm,psfFWHM, fitRingArg{:});
else
    fitParAvg=[];
end

nFr = size(ringStack,3);
ringStack_noBg = 0.*ringStack;
for ii =1:nFr
    
    %fit each ring.
%     display(['Frame: ',num2str(ii)]);
    if ii==1
        [ringIm_noBg, ringKymoCell{ii}, circleData{ii}, rawKymoCell{ii},fitPar(ii,:)] = bgSubAndProfile(ringStack(:,:,ii),pixSzNm,lineWidthNm, psfFWHM,doFixedRadiusFit,fitParAvg,fitRingArg{:});
    else
        [ringIm_noBg, ringKymoCell{ii}, circleData{ii}, rawKymoCell{ii},fitPar(ii,:)] = bgSubAndProfile(ringStack(:,:,ii),pixSzNm,lineWidthNm, psfFWHM,doFixedRadiusFit,fitParAvg,'InitialGuess', fitPar(ii-1,:), fitRingArg{:});
    end
    kymoSz(ii,:) = size(ringKymoCell{ii});
    rNm=circleData{ii}.r*pixSzNm;
    kymoInfo(ii,:) = [ii,rNm,kymoSz(ii,1),kymoSz(ii,2)];
    ringStack_noBg(:,:,ii) = ringIm_noBg;
end

maxKymoWidth = max(kymoSz(:,2));
ringKymograph = zeros(nFr,maxKymoWidth);
rawKymograph = zeros(nFr,maxKymoWidth);
for ii = 1:nFr
    ringKymograph(ii,1:kymoSz(ii,2))=ringKymoCell{ii};
    rawKymograph(ii,1:kymoSz(ii,2))=rawKymoCell{ii};
end

if doZeroPadKymo%no point in padding the raw kymograph
    ringKymograph(end+1,:)=0;
end


%------------------------------------------------------------------
function [ringIm_noBg,ringIntensity, circleData, rawRingIntensity,fitPar] = bgSubAndProfile(ringIm,pixSzNm,lineWidthNm, psfFWHM,doFixedRadiusFit,fitParAvg,varargin)
%extract circular kymograph, integrating over widthNM annulus thickness
%use fitting to the ring to find the diameter, should make it more robust eg on small rings


%fit blurred ring to the image
if doFixedRadiusFit
    [fitPar,~,ringIm_noBg] = fitRing(ringIm, pixSzNm,psfFWHM, 'FixedRadiusFit', fitParAvg, varargin{:});
else 
    [fitPar,~,ringIm_noBg] = fitRing(ringIm, pixSzNm,psfFWHM, varargin{:});
end

x=fitPar(1);
y=fitPar(2);
circ_z=[x,y];
circ_r=fitPar(3);

% calculate a single pixel of circumference and use that as the 
% spacing to sample the circle
theta_step = 1/circ_r;
% non clockwise non-top thetat really confusing
%define clockwise from twelve-o-clock
theta = 0:theta_step:2*pi;
%theta = pi/2:-theta_step:-3/2*pi;
distStep =cumsum([0,ones(1,numel(theta)-1)*abs(theta(2)-theta(1))]);
Pcirc = [circ_z(1)+circ_r*cos(theta(:)), circ_z(2)+circ_r*sin(theta(:)), theta(:), distStep(:)];


%use this to plot a profile for all frames (lineWidth wide)
lineWidthPix = lineWidthNm/pixSzNm;
ringProfile=[];
[ringIntensity] = getProfile(ringIm_noBg,Pcirc,lineWidthPix,circ_z,circ_r,pixSzNm);
ringIntensity = ringIntensity(:)';
[rawRingIntensity] = getProfile(ringIm,Pcirc,lineWidthPix,circ_z,circ_r,pixSzNm);
rawRingIntensity = rawRingIntensity(:)';

circleData.coord = Pcirc;
circleData.z = circ_z;
circleData.r = circ_r;

%------------------------------------------------------------------
function [profileIntensityWide, profileIntensity1pix] = getProfile(I,samplePts,lineWidthPix,circ_z,circ_r,pixSzNm);

x = samplePts(:,1);
y= samplePts(:,2);
nPt = numel(x);

STEPSZ =0.1;
%nLineStep is in each direction +/-
nLineStep = round(lineWidthPix/2/STEPSZ);

%for each point, get the intensity
[Xim,Yim] = meshgrid(1:size(I,2),1:size(I,1));
for ii = 1:nPt
    %get a perpendicular line defining the width to integrate along
    XC = [x(ii),y(ii)];
    dX = [XC(1) - circ_z(1),XC(2) - circ_z(2)];
    dXnorm = dX./sqrt(sum(dX.^2));%normalize
    dL = -nLineStep*.1:.1:nLineStep*.1;
    Xline = [repmat(XC',1,numel(dL)) + dXnorm'*dL]';
    %%interpolate the intensity to sub-pixel precision
    lineProfile = interp2(Xim,Yim,I,Xline(:,1),Xline(:,2),'cubic');
    profileIntensityWide(ii) =trapz(dL,lineProfile);
    profileIntensity1pix(ii) = interp2(Xim,Yim,I,x(ii),y(ii),'cubic');
end


