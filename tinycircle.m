function [plist,particles] = tinycircle(simName)

% initialize the simulation to get its dimensions
sim = MDsim();

% make the shape of the circular subgrain, using sim's dimensions
sim.width = 250;
sim.height = 250;

theta = 0:.1:2*pi;
x = sim.width/2 + 30*cos(theta);
y = sim.height/2 + 30*sin(theta);
sim.subgrain_cell = {[x',y']};
sim.subgrain_orientations = [pi/6];

% initialize the polycrystal
sim = sim.initialize_grains();
sim.preview_setup();

% let polycrystal thermalize for 10 steps TODO fix to 10
sim.checkpoint_frequency = 2;
sim.num_frames  = 10;
[plist, particles] = sim.run_sim(strcat(simName,'/chkpt'));


xlim1 = 10;
xlim2 = 230;
ylim1 = 10;
ylim2 = 230;

% want: csvs corresponding to each frame saved so i can pick the ones i like!
for frameNum = 2:2:10 % TODO stop hardcoding u buffoon
    
    % cut out all the particles within a window
    parts = plist(plist(:,3) == frameNum, :);
    window = plist(plist(:,3) == frameNum & ...
                     plist(:,1) >= xlim1 & ...
                     plist(:,1) <= xlim2 & ...
                     plist(:,2) >= ylim1 & ...
                     plist(:,2) <= ylim2, :);

    newColumn = zeros(size(window,1),1);
    % x, y, a column that will become psi6 but is all 0s for now, and particleID
    data = [window(:,1:2),newColumn,window(:,4)];

    % get |psi6| data
    psi6_frame = psi6_simulation(parts,sim.width,sim.height,sim.particle_diam*sim.tightness);
    for i=1:size(data,1)
        pID = data(i,4);
        % look up this particle's |psi6| by particle ID
        %disp(pID);
        psi6 = psi6_frame(psi6_frame(:,4)==pID,3);
        % record value in 4th column of data
        data(i,3) = psi6;
    end

    % output particle positions to a csv
    fname = strcat(simName,'/',simName,'_',num2str(frameNum),'.csv');
    dlmwrite(fname,[20,20,200,200,5]); % first line: [xmin,ymin,window width, window height, bead radius]
    dlmwrite(fname, data(:,1:3), '-append'); % remaining lines: [x position, y position, |psi6|]
end
end

