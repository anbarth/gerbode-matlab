sim = MDsim();

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

%%%%%%%%%%%%%%%%%% 7.5 %%%%%%%%%%%%%%%%%%%
sim.subgrain_cell = {[x',y']};
sim.subgrain_orientations = [7.5*pi/180];
sim.grain_orientation = -1*7.5*pi/180;

sim = sim.initialize_grains();
sim.preview_setup();

sim.checkpoint_frequency = 2;
sim.num_frames  = 6;
sim.run_sim('readshock4_7p5/chkpt');

%%%%%%%%%%%%%%%%%% 12.5 %%%%%%%%%%%%%%%%%%%
sim2 = MDsim();
sim2.subgrain_cell = {[x',y']};
sim2.subgrain_orientations = [12.5*pi/180];
sim2.grain_orientation = -1*12.5*pi/180;

sim2 = sim2.initialize_grains();
sim2.preview_setup();

sim2.checkpoint_frequency = 2;
sim2.num_frames = 6;
sim2.run_sim('readshock4_12p5/chkpt');

%%%%%%%%%%%%%%%%%% 12.5b %%%%%%%%%%%%%%%%%%%
sim3 = MDsim();
sim3.subgrain_cell = {[x',y']};
sim3.subgrain_orientations = [12.5*pi/180];
sim3.grain_orientation = -1*12.5*pi/180;

sim3 = sim3.initialize_grains();
sim3.preview_setup();

sim3.checkpoint_frequency = 2;
sim3.num_frames = 6;
sim3.run_sim('readshock4_12p5b/chkpt');

%%%%%%%%%%%%%%%%%% 12.5c %%%%%%%%%%%%%%%%%%%
sim4 = MDsim();
sim4.subgrain_cell = {[x',y']};
sim4.subgrain_orientations = [12.5*pi/180];
sim4.grain_orientation = -1*12.5*pi/180;

sim4 = sim4.initialize_grains();
sim4.preview_setup();

sim4.checkpoint_frequency = 2;
sim4.num_frames  = 6;
sim4.run_sim('readshock4_12p5c/chkpt');

%%%%%%%%%%%%%%%%%% 15 %%%%%%%%%%%%%%%%%%%
sim5 = MDsim();
sim5.subgrain_cell = {[x',y']};
sim5.subgrain_orientations = [15*pi/180];
sim5.grain_orientation = -1*15*pi/180;

sim5 = sim5.initialize_grains();
sim5.preview_setup();

sim5.checkpoint_frequency = 2;
sim5.num_frames = 6;
sim5.run_sim('readshock4_15/chkpt');

