% initialize the simulation to get its dimensions
sim = MDsim();
sim.width = 1000;
sim.height = 600;

% make the shape of the circular subgrain, using sim's dimensions
theta = 0:.1:2*pi;
x = sim.width/2 + 250*cos(theta);
y = sim.height/2 + 250*sin(theta);
sim.subgrain_cell = {[x',y']};
sim.subgrain_orientations = [pi/6];

% initialize the polycrystal
sim = sim.initialize_grains();
sim.preview_setup();
%particles = sim.initial_particles;

% output particle positions to a csv
%fname = 'myCircularGrain.csv';
%dlmwrite(fname, particles(:,1:2)); 

sim.checkpoint_frequency = 1;
sim.num_frames  = 100;
[~, particles] = sim.run_sim([simName '/chkpt']);



