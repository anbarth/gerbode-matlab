function particles = flatGBsim(theta,simName)

% initialize the simulation to get its dimensions
sim = MDsim();

% make the shape of the rectangular subgrain, using sim's dimensions
x_top = -1:1:sim.width+1;
y_top = round(2*sim.height/3) * ones(1, length(x_top));

x_bot = sim.width+1:-1:-1;
y_bot = round(sim.height/3) * ones(1, length(x_bot));

y_left = round(2*sim.height/3):-1:round(sim.height/3);
x_left = zeros(1, length(y_left));

y_right = round(sim.height/3):1:round(2*sim.height/3);
x_right = (sim.width+1) * ones(1, length(y_right));

x = [x_top, x_right, x_bot, x_left];
y = [y_top, y_right, y_bot, y_left];

% initialize the polycrystal
sim.subgrain_cell = {[x',y']};
sim.subgrain_orientations = [theta*pi/180];
sim.grain_orientation = -1*theta*pi/180;
sim = sim.initialize_grains();
sim.preview_setup();

% let polycrystal thermalize for 10 steps TODO fix to 10
sim.checkpoint_frequency = 5;
sim.num_frames  = 5;
[~, particles] = sim.run_sim([simName '/chkpt']);

% cut all the particles within a window
xlim1 = 590;
xlim2 = 1210;
ylim1 = 240;
ylim2 = 560;
window = particles( particles(:,1) >= xlim1 & ...
                             particles(:,1) <= xlim2 & ...
                             particles(:,2) >= ylim1 & ...
                             particles(:,2) <= ylim2, :);
                         
newColumn = zeros(size(window,1),1);
% x, y, a column that will become psi6 but is all 0s for now, and particleID
data = [window(:,1:2),newColumn,window(:,7)];

% get |psi6| data
psi6_frame = psi6_simulation(particles,sim.width,sim.height,sim.particle_diam*sim.tightness);
for i=1:size(data,1)
    pID = data(i,4);
    % look up this particle's |psi6| by particle ID
    psi6 = psi6_frame(psi6_frame(:,4)==pID,3);
    % record value in 4th column of data
    data(i,3) = psi6;
end

% output particle positions to a csv
fname = [simName '/' simName '.csv'];
dlmwrite(fname,[600,250,600,300,5]); % first line: [xmin,ymin,window width, window height, bead radius]
dlmwrite(fname, data(:,1:3), '-append'); % remaining lines: [x position, y position, |psi6|]

end

