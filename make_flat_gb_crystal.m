particles = cnt_cell{1}(:,1:2);

simHeight = max(particles(:,2));
simWidth = max(particles(:,1));
y_gb = simHeight/3;

myHeight = 200;
myWidth = 800;
buffer = 20;

% cut all the particles within a window
xlim1 = 0;
xlim2 = myWidth+2*buffer;
ylim1 = y_gb-myHeight/2-buffer;
ylim2 = y_gb+myHeight/2+buffer;
window = particles( particles(:,1) >= xlim1 & ...
                             particles(:,1) <= xlim2 & ...
                             particles(:,2) >= ylim1 & ...
                             particles(:,2) <= ylim2, :);

xlim1_inner = buffer;
xlim2_inner = myWidth+buffer;
ylim1_inner = y_gb-myHeight/2;
ylim2_inner = y_gb+myHeight/2;
                         
make_dots_wraparound(window(:,1),window(:,2),5,"mydots.png",ceil(simWidth),ceil(simHeight));
disp(size(window))
% output particle positions to a csv
fname = "../python/crystals/chrisFlat/30_1.csv";
% first line: radius
dlmwrite(fname,5); 
% second line: window
dlmwrite(fname,[xlim1_inner,ylim1_inner,xlim2_inner,ylim1_inner,xlim2_inner,ylim2_inner,xlim1_inner,ylim2_inner], '-append');
% remaining lines: particle centers
myIDs = ( 1:size(window,1) )';
dlmwrite(fname, [myIDs,window], '-append');
writeNeighbs("chrisFlat/30_1");