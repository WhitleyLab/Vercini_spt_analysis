function [FWHM, intensity, widthResult] = fitProfile(im,width,plotpos_y)
%function [FWHM, intensity, widthResult;] = fitProfile(im,width,plotpos_y)

y = double(im(plotpos_y, :));
h = size(im,2);
x= 1:h;
superGaussian =@(a,x)  a(1) +a(2)*exp(-((x-a(3)).^2 / (2*a(4)^2)).^ a(5) );
initwidth = 5;
initguess = [min(y), max(y), mean(x),initwidth,1];
options = optimoptions('lsqcurvefit','Display','off');
a = lsqcurvefit(superGaussian,initguess,x,y,[],[],options);

xC = a(3);
par_d= a(4);
par_e = a(5);
FWHM = 2*sqrt(2)*par_d*(log(2))^(1/(2*par_e));
x0 = xC-FWHM/2;
x1 = xC + FWHM/2;
widthResult=[xC,x0,x1];
intensity = max(y);
end

