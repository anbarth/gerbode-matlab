im = imread('readshock4_15/chkpt6.png');

figure(20);
imshow(im);
hold on;

xmax = 1915;
ymax = 1216;

% old: too long
%xlim1 = 700;
%xlim2 = 1700;
%ylim1 = 200;
%ylim2 = 600;

% new n improved
xlim1 = 600;
xlim2 = 1200;
ylim1 = 250;
ylim2 = 550;

a1 = [xlim1, xlim1];
b1 = [0, ymax];
a2 = [xlim2, xlim2];
b2 = [0, ymax];
a3 = [0, xmax];
b3 = [ylim1, ylim1];
a4 = [0, xmax];
b4 = [ylim2, ylim2];

plot(a1,b1,'LineWidth',3);
plot(a2,b2,'LineWidth',3);
plot(a3,b3,'LineWidth',3);
plot(a4,b4,'LineWidth',3);