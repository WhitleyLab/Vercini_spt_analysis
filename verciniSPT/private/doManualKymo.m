function [ ringKymograph, circleData] = doManualKymo(ringStack,pixSzNm)

ringMax = max(ringStack,[],3);
h=figure;
imagesc(ringMax);
colormap('gray');
axis equal;
display('draw a circle and then press any key');
display('you can drag and resize if you need to');
roi=drawcircle();
pause;
display('calculating kymograph');
circ_z= roi.Center;
circ_r = roi.Radius;
close(h);

lineWidthNm = pixSzNm;
nFr = size(ringStack,3);
for ii =1:nFr
    ringIm  = ringStack(:,:,ii);
    [ringKymograph(ii,:), circleData] = getKymo(ringIm,pixSzNm,lineWidthNm,circ_z,circ_r);
end



%---------------------------------------------------------------------------------
function [ringIntensity, circleData] =getKymo(ringIm,pixSzNm,lineWidthNm, circ_z,circ_r)
%x=fitPar(1);
%y=fitPar(2);
%circ_r=fitPar(3);

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
[ringIntensity] = getProfile(ringIm,Pcirc,lineWidthPix,circ_z,circ_r,pixSzNm);
ringIntensity = ringIntensity(:)';

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
%     I = double(I); %% HACK!!!!!!!! 220831 kw
    lineProfile = interp2(Xim,Yim,I,Xline(:,1),Xline(:,2),'cubic');
    profileIntensityWide(ii) =trapz(dL,lineProfile);
    profileIntensity1pix(ii) = interp2(Xim,Yim,I,x(ii),y(ii),'cubic');
end


