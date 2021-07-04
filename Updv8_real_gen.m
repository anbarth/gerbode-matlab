function [ matrix_data,right_wall,top_wall ] = Updv8_real_gen( matrix_data, particle_size, neighb_threshold, step_size, frame_number, incoming_rays,right_wall,top_wall, laser_on, laser_off, force_threshold, laser_origin)
%This program is a variation of the Upd8, except time units are scaled by
%seconds instead of being unitless.  This is done by using a Brownian
%dynamics approach to simulation.  Wiki it if more details
%needed.Essentially we use the diffusion constant to approximate a velocity
%at each time step and integrate normally.
%The inputs are:
%matrix_data: matrix containing all the information of our particles.  nx7
%matrix with  columns x-y position, x-y veloity, x-y acceleration, index
    %number
%particle size: diameter of our particles in pixels
%neighb_threshold: the number of particle diameters away to check for
%neighboring particles
%step_size: the time step in seconds
%frame_number: the current frame number from higher level programs
%incoming_rays: a matrix detailing the origins and direction light rays
%   which refract through our colloids to impart a force
%right_wall: the x-coordinate of our boundary condition, with the
%left_wall: zero
%top_wall: the y coordinate of our boundary conditioin with bottom wall
%   being zero
%laser_on: the frame at which the laser turns on
%laser_off: the frame at which the laser turns off
%force_threshold: the number of particle diameters away to calculate laser
%   forces for (if this is too big it will Not End Ever)
%laser_origin: the location of the laser (optional), typically [1160 450]

%Jeremy Wang '17, edited by Nina Brown '19



%Define some constants
%mass=1.38*10^-15;
viscosity=1.99*10^-3;
boltz_const = 1.381*(10^(-23));
temp=300;
%scaling_factor is pixels/meter. We want each particle to have a diameter
%of 1.2 microns, or 1.2 * 10^(-6) m
scaling_factor=particle_size/(1.2*(10^-6));
%refractive indices for DMSO and particles, respectively
n1 = 1.477; 
n2 = 1.424;
%where's the laser? it's an optional input now
if nargin < 12
    laser_origin = [1160 450];
end

%determine when the laser will be on
%Note that if laser_on = laser_off, the trap is never on.
if frame_number >= laser_on && frame_number < laser_off
    trap='on';
else
    trap='off';
end


%create neighbors list that has the neighboring particles for each particle
neighbors_list=Find_neighbors_wraparound(matrix_data,particle_size, neighb_threshold,right_wall,top_wall);


%repeat trajectory calculations for a specified amount of times with the
%same neighbor list
for i= 1:10
    
    %From Stokes drag
    %squig = 6 pi eta r
    %r = particle_size/2/scaling_factor
    %   = particle_size/2/(particle_size/(1.2*(10^-6)))
    %   = 0.6*10^-6 (1.2 micron diameter particles)
    squig=6*pi*viscosity*(particle_size/2/scaling_factor);
    %calculate laser force
    if strcmp(trap,'on')
        %matrix of laser force for each particle
        laser_force=zeros(length(matrix_data(:,1)),2);

        %the laser is centered at coordinates [1160 450], so shift all the
        %particles by the appropriate amount
        sphere_trans=matrix_data(:,[1 2])- repmat(laser_origin, length(matrix_data(:,1)),1);
        %find the distance of the particle from the center of the laser
        sphere_dist=(sphere_trans(:,1).^2 + sphere_trans(:,2).^2).^(1/2);
        
        %find all particles ("colloids" [sic]) within the effect of the laser
        colloids_in_laser=matrix_data(sphere_dist<force_threshold*particle_size,[1 2 7]);
        
        %create matrix with only particles being affected, with their
        %distance from the laser center in pixels as the [x y] at [1 2]
        spheres_matrix=[colloids_in_laser(:,[1 2]) zeros(length(colloids_in_laser(:,1)),1) repmat(particle_size,length(colloids_in_laser(:,1)),1)]...
            -repmat([laser_origin 0 0],length(colloids_in_laser(:,1)),1);
        
        %put in units of microns
        spheres_matrix=spheres_matrix/(scaling_factor*(10^(-6)));
        
        %get the force of the laser's rays on the particles
        [x,y]=Calculate_Momentum(spheres_matrix,n1,n2,incoming_rays);
        
        laser_force(colloids_in_laser(:,3),:)= [x y]                                                                                                                                                                                                                                ;
    else
        %if trap isn't on, the laser forces are zero!
        laser_force=zeros(length(matrix_data(:,1)),2);
    end

    
    %update trajectory using brownian dynamics equation of motion (here 1D):
    % x(t + dt) = x(t) = F_laser / (6 pi eta R) * dt + dx_brownian
    % TODO 1/3/2020 fix the diffusion constant? like it should probably be 2D?
    matrix_data(:,[1 2]) = matrix_data(:, [1 2]) + laser_force(:,[1 2])/squig*scaling_factor*step_size + normrnd(0,sqrt(2*(boltz_const*temp)/squig*step_size),[length(matrix_data(:,1)) 2])*scaling_factor;
    
    %resolve any collisions
    matrix_data= Apply_collision_wraparound(matrix_data, particle_size,neighbors_list, right_wall,top_wall);
    
    %random_force=diffusion/step_size*normrnd(0,1,[length(matrix_data(:,5)) 2]);
    
    %enforce wraparound boundary conditions
    left_border=matrix_data(:,1)<0;
    matrix_data(left_border,1)=matrix_data(left_border,1)+right_wall;
    
    right_border=matrix_data(:,1)>right_wall;
    matrix_data(right_border,1)=matrix_data(right_border,1)-right_wall;
    
    bottom_border=matrix_data(:,2)<0;
    matrix_data(bottom_border,2)=matrix_data(bottom_border,2)+top_wall;
    
    top_border=matrix_data(:,2)>top_wall;
    matrix_data(top_border,2)=matrix_data(top_border,2)-top_wall;
    

end
end

