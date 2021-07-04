function [particles,plist] = No_Draw_Simulate_gen(num_particle,particle_size, neighb_threshold, step_size, right_wall, top_wall, matrix_data, end_num, num_frames, laser_on, laser_off, laser_origin, saving_frequency)
%See No_Draw_Simulate for more detail. 
%-- Nina Brown 2017 but almost all of this code is copied from Jeremy. This
%version takes in number of frames to simulate, and which laser frames
%should be one, so that we don't have to change it within the function.
%INPUTS (copied from N_D_S):
%num particle is the number of particles wanted in our simulation

%particle_size is the diameter of our particles in pixels

%diffusion is the strength of our random Brownian force.  This is unused
%currently, as it is calculated from theory

%friction is the coefficient of friction for the drag force. This is also
%unused as it is also calculated from theory

%gravity is the coefficient for an effectly gravity force.  This points in
%the x direction.  Unused

%threshold is a multiplier that specifies how many particle diameters out we
%wish to check for potential neighbors for each particle

%stepsize is the amount of time each step or frame of the simulation will be

%right wall is the location of the right boundary where wraparound
%conditions will be applied

%top wall is the location of the upper boundary where wraparound conditions
%will be applied

%matrix_data is a preinitialized set of particle coordinates, usually used
%when initilizing grain boundaries

%end_num is the number that will be added at the end of saved frames of the
%simulation

%num_frames: the number of steps/frames to run the simulation for

%laser_on: the frame at which the laser should turn on

%laser_off: the frame at which the laser should turn off

%laser_origin is the [x y] position of the laser

%saving_frequency is how frequently (every x number of frames)
%we want to append our particles to
%plist. For long simulations this will prevent the size of plist from being
%too large.

%OUTPUTS:
%particles: list of all particles in each frame
%plist: list of ID'd particles in each frame
%gb_cell: cell containing list of gb location for each frame

%Jeremy Wang '17, edited by Nina Brown '19

%optional inputs
if nargin < 13
    saving_frequency = 1;
end   
if nargin < 12
    laser_origin = [1160, 450];
end   

%initialize a counter
counter=1;

%specify how many steps the simulation will run for
steps=num_frames;

%initialize a plist, which will hold particle positions for each frame of
%the simulation, and gb_cell, which will hold the GB position for each
%frame
frames_saved = floor(steps/saving_frequency); 
plist = zeros(num_particle*frames_saved,4);

%initialize rays that will be used to calculate laser force
%max radius from laser to calculate rays = 9 microns
%d = 0.15 microns (see JW thesis) (somewhere b/w 0.1-0.9??)
%total power = 0.7 W (~max power of our laser)
%beam waist = 4????? 0.13????? (in JW code this is beam_spot but it's actually the
%   beam waist w(0) instead of w(z)) (we calculated this to be 0.13 but
%   this seems small...) so i guess 4 sure
%   Update: estimated w/ hole diameter from a very loose ramp cell: 4.5 microns
%   UPDATED AGAIN 15 microns WhOoPs!
%   Okay so 15 microns makes this essentially unusable and also would
%   probably make holes pretty big (? idrk i tried to run this and it
%   didn't finish in like 5 days so changing it back to like 4 bc that
%   seems to give reasonably sized holes
ray_radius = 9;
d_from_focus = 0.15;
laser_power = 0.7;
beam_waist = 4;
if laser_on ~= laser_off
    incoming_rays=circular_ray_make(ray_radius, d_from_focus, laser_power, beam_waist);
else
    incoming_rays = zeros(10,7);
end
%set the radius for which we calculate forces on the particles to be the
%same as the radius for which we calculated rays
force_threshold = ray_radius;

%check how many particles are in the simulation currently
current_num=length(matrix_data(matrix_data(:,1)>0 & matrix_data(:,2)>0)) + 1;

%simulation loop
for i=1:steps
    disp_string = strcat('Starting frame ', num2str(i), ' in Simulation ', num2str(end_num));
    %add particles if we do not have the desired number of particles
    if current_num<num_particle
        x_co=right_wall*rand;
        y_co=top_wall*rand;
        while sum(matrix_data(:,1)<x_co + 1.2*particle_size &  matrix_data(:,1)>x_co - 1.2*particle_size ...
                & matrix_data(:,2)<y_co + 1.2*particle_size &  matrix_data(:,2)>y_co - 1.2*particle_size)~=0
            x_co=right_wall*rand;
            y_co=top_wall*rand;
        end
        matrix_data(current_num,[1 2])=[x_co y_co];
        current_num=current_num + 1;
    end
    
    %make a name 
    name=strcat('sim_',num2str(end_num),'_frame_', num2str(i),'.png');
    
    %here we check how many of our particles are not 0 entries.  In
    %essence how many actual particles we currently have on screen
    
    non_zero=matrix_data(:,1)>0;
    
    %gets rid of zero entries in our matrix
    non_zero_matrix=matrix_data(matrix_data(:,1) ~=0 & matrix_data(:,2) ~= 0,:);

    %{
    %draw a frame with particles using make_dots_wraparound
    make_dots_wraparound(non_zero_matrix(:,1),non_zero_matrix(:,2),particle_size/2,name,right_wall,top_wall);
    %}
    string=strcat('frame ', num2str(i),' done');
    disp(string)
    
    
    %update particle positions using Updv8_real_gen
    [non_zero_matrix,right_wall,top_wall] = Updv8_real_gen( matrix_data, particle_size, neighb_threshold, step_size, i, incoming_rays, right_wall, top_wall, laser_on, laser_off, force_threshold, laser_origin);
    

    %add zeros back in for consistency
    matrix_data(1:length(non_zero_matrix(:,1)),:)=non_zero_matrix;
    
    %Don't always append to plist
    if mod(i,saving_frequency) == 0
        plist(counter:counter+num_particle -1, 1)=matrix_data(:,1);
        plist(counter:counter+num_particle -1, 2)=matrix_data(:,2);
        plist(counter:counter+num_particle -1,3)=i*ones(num_particle,1);
        plist(counter:counter+num_particle -1,4)=[1:num_particle]';
        %update counter
        counter=counter+num_particle;
    end
    
end
%output matrix_data
particles = matrix_data;

%these commands save particles and plist if desired
%filename=strcat( 'Particle_Positions_Colloid_',num2str(diffusion), '_', num2str(friction), '_', num2str(gravity),'_',end_str,'.mat');
%save(filename, 'particles');

%filename2=strcat( 'Particle_Plist_Colloid_',num2str(diffusion), '_', num2str(friction), '_', num2str(gravity),'_',end_str,'.mat');
%save(filename2, 'plist');
end
