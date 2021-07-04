xlim1 = 790;
xlim2 = 1210;
ylim1 = 240;
ylim2 = 560;

bnt = particles( particles(:,1) >= xlim1 & ...
                             particles(:,1) <= xlim2 & ...
                             particles(:,2) >= ylim1 & ...
                             particles(:,2) <= ylim2, :);

                         
simName = 'readshock4_15';
%writematrix(bnt, [simName '/' simName '.csv']);
dlmwrite('csvTest.csv',[400,300,5.5]);
dlmwrite('csvTest.csv', bnt, '-append');