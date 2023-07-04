function plotCircleFig(im, circleData)

imagesc(im);
colormap gray;
hold all;
%fit
x=[circleData.coord(:,1);circleData.coord(1,1)];
y=[circleData.coord(:,2);circleData.coord(1,2)];
plot(x,y,'r--','LineWidth',1.5);
axis equal;
