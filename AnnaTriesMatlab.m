im = imread('200222A_ts_t026.png');
d_im = double(im);
b = bpass(d_im, 1, 20);
pk = pkfnd(b,1,20);

%figure(26);
%imshow(b);

figure(27);
imshow(im);
hold on;
x = pk(:,1);
y = pk(:,2);
%scatter(x,y,'b','filled');

% original myParts that i did lots of stuff w
% x bounds: 550 to 900
% y bounds: 400 to 700
% a more orderly region
% x bounds: 150 to 400
% y bounds: 400 to 700
xlim1 = 150;
xlim2 = 500;
ylim1 = 400;
ylim2 = 700;

a1 = [xlim1, xlim1];
b1 = [0, 1500];
a2 = [xlim2, xlim2];
b2 = [0, 1500];
a3 = [0, 2000];
b3 = [ylim1, ylim1];
a4 = [0, 2000];
b4 = [ylim2, ylim2];


plot(a1,b1,'LineWidth',3);
plot(a2,b2,'LineWidth',3);
plot(a3,b3,'LineWidth',3);
plot(a4,b4,'LineWidth',3);

cnt = cntrd(b,pk,25);
bnt = cnt( cnt(:,1) >= xlim1 & cnt(:,1) <= xlim2 & cnt(:,2) >= ylim1 & cnt(:,2) <= ylim2, :);
 
writematrix(bnt, 'myPartsOrderly.csv');
success = fclose('all');



