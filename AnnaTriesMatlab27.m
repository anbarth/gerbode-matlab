im = imread('200222A_ts_t027.png');
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


% a more orderly region
% bnt is 213
%xlim1 = 150;
%xlim2 = 500;
%ylim1 = 300;
%ylim2 = 600;

% a more orderly region
% volfrac 0.376
%xlim1 = 150;
%xlim2 = 600;
%ylim1 = 800;
%ylim2 = 1200;

% a less orderly region
% volfrac 0.388
xlim1 = 1300;
xlim2 = 1700;
ylim1 = 500;
ylim2 = 900;

% a less orderly region
% bnt is 209
%xlim1 = 750;
%xlim2 = 1100;
%ylim1 = 200;
%ylim2 = 500;

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


%writematrix(bnt, 'myPartsOrderly.csv');
%success = fclose('all');



