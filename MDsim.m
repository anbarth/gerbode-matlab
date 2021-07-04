classdef MDsim
    properties
      %-------------
      % Simulation parameters
      %-------------
        particle_diam %particle diameter
        height %height of our simulation window
        width %width of our simulation window
        neighb_threshold %distance where we check for neighbors to detect collisions 
        step_size %simulation timestep in seconds
        num_frames %number of frames we are simulating. Each frame is 10 timesteps
        saving_frequency %how frequently do we want to save data in plist
        
        %an nx7 list that serves as the initial conditions
        %for the simulation.
        %Each row is
        %[x y x_vel y_vel x_accel y_accel particleID]
        %Velocity and accceleration are all initialized as 0.
        initial_particles
                          
      
      %-------------
      % Stuff for initializing the first frame
      %-------------
      
        %The center to center distance, as a function of particle diameter
        tightness
      
        %The orientation
        %(defined as smallest angle between a neighbor and the horozontal)
        %of the 'main' grain
        grain_orientation
      
        %an nx1 cell where each entry is an nx2 array representing 
        %vertices(x y) of a polygon that completely bounds a subgrain
        subgrain_cell 
      
      
        %an nx1 array where each entry represents the orientation
        %(defined as smallest angle between a neighbor and the horozontal)
        %of the subgrains
        subgrain_orientations 

      %-------------
      % Stuff for resuming simulations
      %-------------
      
        %the current frame number
        current_frame
        
        %the current frame's particles
        current_particles
        
        %saves a plist and particles after simulating this many frames.
        %This may be changed slightly so that it is a multiple of
        %saving_frequency.
        checkpoint_frequency
      
        
        

      %-------------
      % Laser Properties
      %-------------      
        laser_on %an integer that says which frame the laser turns on
        laser_off %an integer that says which frame the laser turns off
        laser_origin %a 1x2 vector (x y) representing the origin of the laser
        
    end
    
    methods
   % ==================================================
   % Constructors
   % ==================================================
       %Default constructor makes a 350-frame long simulation the a box 
       %the same size as our real experiment images,
       %with reasonable choices of parameters.
        function self = MDsim()
            if nargin == 0
                self.particle_diam = 10; 
                self.height = 1216;
                self.width = 1920;
                self.neighb_threshold = 2; 
                self.step_size= .005;
                self.num_frames = 350;
                self.saving_frequency = 1;
                
                self.tightness = 1.01; 
                self.grain_orientation = 0;
                self.subgrain_cell = {};
                self.subgrain_orientations = [];
                
                self.current_frame = 0;
                self.checkpoint_frequency = 250;

                self.laser_on = 0;
                self.laser_off = 0;
                self.laser_origin = [self.width/2, self.height/2];         
            else
                error('Incorrect number of inputs for constructor')
            end
        end
      
   % ==================================================
   % Grain Initializing
   % ==================================================   
        %Create the specified grains and creates
        %an nx7 list that serves as the initial conditions
        %for the simulation.
        %Each row is [x y xvel yvel xacc yacc particleID]
        %
        %Note that this returns a new MDsim object. This new object is the
        %one that has the initialized grain.
        %
        %Inputs:
        %addParticles: boolean if the function addParticles should be run
        %to incease the particle density of the inner grain - default false
        %numParticles: the number of particles to add - default 0
        %buffer: particles will add in places with 'buffer'*particle diameter amount of space
        %between particles - default 0.9
        %optimize_placements gives the option to try slightly different
        %placements of the subgrains to maximize the number of particles in
        %the simulation
        function self = initialize_grains(self, addParticles, numParticlesToAdd, buffer, optimize_placements)
            
            if nargin < 2
                addParticles = false;  %%NEWLY ADDED TO ACCOUNT FOR ADDING PARTICLES IN SIMS RESUME SCRIPTS
            end
            
            if nargin < 3
                numParticlesToAdd = 0;
            end
            
            if nargin < 4
                buffer = 0.9;
            end
            
            if nargin < 5
                optimize_placements = false;
            end
            
           
            %Get the main grain
            [self, particles, pSize, width_] = self.makeMainGrain();
            disp(strcat("Particle size is now ",num2str(self.particle_diam)));
            disp(strcat("Width is now ", num2str(self.width)));
            self.particle_diam = pSize;
            self.width = width_;
            %best_particles = particles; % TODO i added this line
            for i = 1:length(self.subgrain_cell)
                %Make a lattice at the orientation of the subgrain
                sub_grain = make_hcp(self.tightness*self.particle_diam/2,...
                                    self.height/2,...
                                    self.width/2,...
                                    ceil(1.2*self.width/self.particle_diam),...
                                    self.subgrain_orientations(i));
                
                %If we are optimizing the placements, come up with a small
                %displacement to the center of the grain
                x_disp_to_try = 0;
                y_disp_to_try = 0;
                if optimize_placements
                    x_disp_to_try = round(self.particle_diam/2);
                    y_disp_to_try = round(self.particle_diam/2);
                end
                
                best_disp = [0 0];
                best_n = 0;
                best_particles = particles;
                for x_disp = -x_disp_to_try:0.25:x_disp_to_try
                    for y_disp = -y_disp_to_try:0.25:y_disp_to_try
                        displacement = [x_disp y_disp];
                        %displace the boundary
                        boundary = self.subgrain_cell{i} + displacement;
                        %Keep the main grain particles outside of the polygon
                        toKeep = ~inpolygon(particles(:,1), particles(:,2), boundary(:,1), boundary(:,2));
                        %Place subgrain particles inside of the poygon
                        toPlace = inpolygon(sub_grain(:,1), sub_grain(:,2), boundary(:,1), boundary(:,2));

                        %Now combine the two.
                        %There are going to be some overlapping particles which we
                        %will need to resolve
                        test_particles = [sub_grain(toPlace,:);particles(toKeep,:)];
                        labels = zeros(1,length(test_particles));
                        labels(1:length(sub_grain(toPlace,:))) = 1;
                        [test_particles, ~] = self.resolveCollisions(test_particles,labels);

                        %There also may be some overlapping particles from the
                        %boundary conditions
                        [test_particles, ~] = self.resolveBoundaryCollisions(test_particles);
                        
                       
                        
                        nParticles = size(test_particles, 1);
                        %check to see if this run was the best
                        if nParticles > best_n
                            best_disp = displacement
                            best_particles = test_particles;
                            best_n = nParticles
                        end
                    end
                end
            
                particles = best_particles; % TODO i added this
            end
            
            %particles = best_particles; % TODO see 12/18
            
            %Now fill in the rest of the columns
            particles = [particles zeros(length(particles(:,1)),4)];
            particles = [particles (1:length(particles(:,1)))'];
            self.initial_particles = particles;
            self.current_particles = particles;
            
            %Now add the particles
            if addParticles == true
                newParticles = self.add_particles(numParticlesToAdd, buffer);
                self.initial_particles = newParticles;
                self.current_particles = newParticles;
            end
            
        end
        
        
        %Makes the background grain. Tries a host of close by particle
        %sizes and widths to get the best fitting possible boundary
        %conditions. Note that in the case of the background orientation
        %being 0, I made special calculations so it can work well.
        %                                                             
        %INPUTS:
        %    particles: an nx2 matrix of particles positions [x y]
        %
        %OUTPUTS:
        %    particles: an nx2 matrix of particle positions [x y] filtered
        %    so that 
        %
        
        
        function [self, mainGrain, newPSize, newWidth] = makeMainGrain(self)
            
            %If it's not the special 0 case
            %if self.grain_orientation ~= 0
            if true % hahaha sorry its me again TODO 1/17
                %make a bunch of nearby sizes and widths to see which one is best
                close_sizes = linspace(.98*self.particle_diam, 1.02*self.particle_diam, 20);
                close_widths = linspace(self.width-round(self.particle_diam)/2,self.width+round(self.particle_diam)/2,round(self.particle_diam) + 1);
                
                maxVolFrac = 0;
                bestW = -1;
                bestD = -1;
                
                for sizeInd = 1:length(close_sizes)
                    for widthInd = 1:length(close_widths)
                    
                    test_diam = 10; % TODO back on my bullshit 1/17
                    %test_diam = self.particle_diam; % TODO this is anna's naughtiness 1/16
                    %test_diam = close_sizes(sizeInd);
                    test_width = close_widths(widthInd);
                    
                    %makes a hcp lattice to match the 'main' grain. Make
                    %sure it's still building centered on the same spot
                    hcp = make_hcp(self.tightness*test_diam/2,...
                                  self.height/2,self.width/2, ...
                                  ceil(1.2*test_width/test_diam),...
                                  self.grain_orientation);

                    %keep particles that are in the image
                    in_image = (hcp(:,1)> 0 & ...
                               hcp(:,1) < test_width & ...
                               hcp(:,2)> 0 & ...
                               hcp(:,2) < self.height );
                    particles = hcp(in_image,:);
                   
                    %resolve collisions due to periodic bc's and check how
                    %good of a fit it was
                    oldWidth = self.width;
                    oldDiam = self.particle_diam;
                    self.width = test_width;
                    [particles, ~] = self.resolveBoundaryCollisions(particles);
                    volFrac = self.getVolFrac(particles);
                    self.particle_diam = oldDiam;
                    self.width = oldWidth;
                    %save the results if they are good
                    if volFrac > maxVolFrac
                        maxVolFrac = volFrac;
                        bestW = test_width;
                        bestD = test_diam;
                    end
                    
                    
                    end
                end
                   
                %Save the mins and recreate the lattice
                newPSize = bestD;
                newWidth = bestW;
                
                %makes a hcp lattice to match the 'main' grain
                hcp = make_hcp(self.tightness*test_diam/2,...
                              self.height/2,newWidth/2, ...
                              ceil(1.2*newWidth/test_diam),...
                              self.grain_orientation);

                %keep particles that are in the image
                in_image = (hcp(:,1)> 0 & ...
                           hcp(:,1) < newWidth & ...
                           hcp(:,2)> 0 & ...
                           hcp(:,2) < self.height );
                particles = hcp(in_image,:);

                %resolve collisions due to periodic bc's
                self.width = newWidth;
                self.particle_diam = newPSize;
                [mainGrain, ~] = self.resolveBoundaryCollisions(particles);
                
                
        
            else
                
                %What size particle so that the top row has no overlap
                repsTall = floor(self.height/(self.tightness*sqrt(3)*self.particle_diam));
                newPSize = self.height/(self.tightness*sqrt(3)*repsTall);
                
                %optimize the width to get the minimum overlap
                %horizontally. Try widths 5 on each side, and get the one
                %with the highest volfrac
                maxVolFrac = 0;
                bestW = -1;
               

                close_widths = linspace(self.width-5,self.width+5,11);
                
                for widthInd = 1:length(close_widths)

                    test_width = close_widths(widthInd);

                    %makes a hcp lattice to match the 'main' grain
                    hcp = make_hcp(self.tightness*newPSize/2,...
                                  self.height/2,self.width/2, ...
                                  ceil(1.2*self.width/newPSize),...
                                  self.grain_orientation);

                    %keep particles that are in the image
                    in_image = (hcp(:,1)> 0 & ...
                               hcp(:,1) < test_width & ...
                               hcp(:,2)> 0 & ...
                               hcp(:,2) < self.height );
                    particles = hcp(in_image,:);

                    %resolve collisions due to periodic bc's and check how
                    %good of a fit it was
                    oldWidth = self.width;
                    oldDiam = self.particle_diam;
                    self.width = test_width;
                    [particles, ~] = self.resolveBoundaryCollisions(particles);
                    volFrac = self.getVolFrac(particles);
                    self.particle_diam = oldDiam;
                    self.width = oldWidth;
                    %save the results if they are good
                    if volFrac > maxVolFrac
                        maxVolFrac = volFrac;
                        bestW = test_width;
                    end


                end
                
                %get the best one out
                newWidth = bestW;
                
               
                %makes a hcp lattice to match the 'main' grain. Make sure
                %to keep it centered at the old place
                hcp = make_hcp(self.tightness*newPSize/2,...
                              self.height/2,self.width/2, ...
                              ceil(1.2*self.width/newPSize),...
                              self.grain_orientation);

                %keep particles that are still in the image
                in_image = (hcp(:,1)> 0 & ...
                           hcp(:,1) <= newWidth & ...
                           hcp(:,2)> 0 & ...
                           hcp(:,2) <= self.height);
      
                 mainGrain = hcp(in_image,:);
                 
                 %Get rid of overlaps
                 oldWidth = self.width;
                 self.width = newWidth;
                 self.particle_diam = newPSize;
                 [mainGrain, ~] = self.resolveBoundaryCollisions(mainGrain);
                 disp('The size of mainGrain after makeMainGrain function is... ');
                 disp(size(mainGrain));


            end
            
                   
            
        end
        
         %Resolves collisions due to the wraparound boundary conditions
        %                                                             
        %INPUTS:
        %    particles: an nx2 matrix of particles positions [x y]. Assumes
        %    that all of them are centered within the frame
        %
        %OUTPUTS:
        %    particles: an nx2 matrix of particle positions [x y] filtered
        %    so that no particles overlap due to the periodic boundary
        %    conditions
        %
        %    removed: the number of particles that were removed in the
        %    process
        %
        
        function [particles, pRemoved] = resolveBoundaryCollisions(self, particles)
            %keep count of particles removed
            pRemoved = 0;
            
            %now get those that are close to the top boundary
            nearTop = (particles(:,2)> self.height - 2*self.particle_diam);
            
            %Move the particles away from the wraparound, resolve
            %collisions and then move them back
            shifted = particles;
            shifted(:,2) = mod(particles(:,2) + self.height/2,self.height);
            [particles, removed1] = self.resolveCollisions(shifted, nearTop);
            pRemoved = pRemoved + removed1;
            particles(:,2) = mod(particles(:,2) - self.height/2,self.height);
            
            %repeat for the right boundary
            nearRight = (particles(:,1)> self.width - 2*self.particle_diam);
            
            shifted = particles;
            shifted(:,1) = mod(particles(:,1) + self.width/2,self.width);
            [particles, removed2] = self.resolveCollisions(shifted, nearRight);
            pRemoved = pRemoved + removed2;
            particles(:,1) = mod(particles(:,1) - self.width/2,self.width);
            
            
        end
        
        
        %Resolves collisions between particles in two different sectors,
        %not taking periodic boundary conditions into account
        %                                                             
        %INPUTS:
        %    particles: an nx2 matrix of particles positions [x y] 
        %    labels: an nx1 matrix that is either 0 or 1, such that no 0
        %    particles overlap with any 1 particles.
        %
        %OUTPUTS:
        %    particles: an nx2 matrix of particle positions [x y] filtered
        %    so that none of the particles overlap in the body (i.e. not
        %    due to the periodic boundary conditions).
        %
        %    totalRemoved: the number of particles that were removed in the
        %    process
        %
        
        function [particles, totalRemoved] = resolveCollisions(self, particles, labels)
            
            %Go through all the particles in the smaller grain and remove
            %all nearest neighbors that are both from the other grain
            %and overlap.
            %Keep repeating the process until there is no overlap
            nParticles = length(labels);
            n0 = nParticles-sum(labels);
            n1 = sum(labels);
            %Let 1 be the smaller grain, so swap labels if 1 is bigger
            if n1 > n0
                labels = ~labels;
                n = n1;
                n1 = n0;
                n0 = n;
            end
            
            pIndices = find(labels);
            removed = 777;
            totalRemoved = 0; %keeps track of the total # of particles removed
            newParticles = particles;
            while removed ~= 0
                %Do a Triangulation to get NNS
                TRI = delaunayTriangulation(newParticles(:,1), newParticles(:,2));
                %initialize the stuff that we'll remove
                toRemove = zeros(1,n0);
                removed = 0;
                %Go through all particles in grain 1
                for pNum = 1:n1
                     % Grab the particle of interest and its nearest neighbors:
                     pIndex = pIndices(pNum);
                     particle = newParticles(pIndex,:); 
                     nns = nearest_neighbors(TRI, pIndex);
                     %now go through its nearest neighbors and see if we are too close
                     for nNum = 1:length(nns)
                         %If the particle is in grain 0 check for overlaps
                         if ~labels(nns(nNum))
                             neighbor = newParticles(nns(nNum),:);
                             distance = norm(particle-neighbor);
                             if distance < self.particle_diam 
                                 toRemove(removed+1) = nns(nNum);
                                 removed = removed +1;
                             end
                         end
                     end
                end

                %Do some logical indexing to remove the marked
                %particles
                if removed > 0
                    toRemove = toRemove(1:removed);
                    removeIndices = zeros(1,length(newParticles));
                    removeIndices(toRemove) = 1;
                else
                    removeIndices = zeros(1,length(newParticles));
                end

                keepIndices = ~removeIndices;
                newParticles = newParticles(keepIndices,:);
                labels = labels(keepIndices);
                pIndices = find(labels);
                
                %and add to the count of removed particles
                totalRemoved = totalRemoved+removed;

            end

            %update our particle list when we've removed all that we
            %need to
            particles = newParticles;
        end
        
        
        %add_particles adds in particles to the region vacated by the two
        %grains, allowing us to even out the total number of particles in
        %simulations with different grain shapes.
        %
        %
        function newInitialParticles  = add_particles(self, numToAdd, buffer, stepSize)
            
            if nargin<2
                numToAdd = 0;
            end
            
            if nargin<3
                buffer = 0.9;
            end
            if nargin<4
                stepSize = 20;
            end
            
            numAdded = 0;
            % First define the list of all particles in the simulation
            allParticles = [self.initial_particles(:,1),self.initial_particles(:,2)];
            
            % Go through the subgrains (embedded grains) and define grids
            % of positions near the grain where we could add particles
            for subgrain = 1:length(self.subgrain_cell)
                
                % Get list of particle positions in the subgrain
                subGB = self.subgrain_cell{subgrain};
                
                % Get the edges of the region where those particles are
                % located
                right = max(subGB(:,1));
                left =  min(subGB(:,1));
                top = max(subGB(:,2));
                bottom = min(subGB(:,2));
                
                % Define a list of potential x and y positions where we could add
                % particles in that region
                pot_x = left-self.particle_diam:self.particle_diam/stepSize:right+self.particle_diam;
                pot_y = bottom-self.particle_diam:self.particle_diam/stepSize:top+self.particle_diam;
                                
                % Loop through a grid defined by pot_x and pot_y and try to
                % find spaces to add a particle
                for ix = 1:length(pot_x)
                    for iy = 1:length(pot_y)
                        
                        % Find the list of distances from the current grid
                        % position to the positions of all particles
                        distances = sqrt((allParticles(:,1) - pot_x(ix)).^2 + (allParticles(:,2) - pot_y(iy)).^2);
                        
                        % Find the distance to the nearest particle
                        dmin = min(distances);
                        if (dmin > self.particle_diam*buffer) && (numAdded < numToAdd) %% Is this where it should be???
                            numAdded = numAdded + 1;
                            % Update self.initial_particles to add in a new
                            % particle at that grid position.
                            self.initial_particles = [self.initial_particles; pot_x(ix) pot_y(iy) 0 0 0 0 (length(self.initial_particles)+1)];
                            % Update allParticles for the next loop
                            allParticles = [self.initial_particles(:,1),self.initial_particles(:,2)];
                        end
                    end
                    
                end 
                
                
                
                
                
            end
            disp(['Yay! ' num2str(numAdded) ' particles added!']);
            newInitialParticles = self.initial_particles;
        end
        
        
        %Thermalize takes self.initial_particles and runs the first
        %timestep of the simulation in a series of much smaller steps. This
        %is equivalent to running the simulation at a lower temperature for
        %the first steps and slowly ramping it up.
        %
        %Inputs:
        %   nSteps: the number of steps to take.
        %   expo: exponent for the ramp function
        %
        %
        function particles = thermalize(self,nSteps, expo)
            %played around with a couple different scalings. Seems like
            %we want one that makes the total thermalizing time less.
            %^4 was a sorta arbitrary choice
            if nargin < 3
                expo = 4;
            end
        
            %make a list of the timesteps to take
            dt_list = zeros(1,nSteps);
            for step = 1:nSteps

                dt_list(step) = (step/nSteps)^expo;
            end
            %now scale it by the actual step size
            dt_list = dt_list.*self.step_size;
            
            %Run simulation for the specified timesteps
            old_step_size = self.step_size;
            old_num_frames = self.num_frames;
            old_saving_frequency = self.saving_frequency;
            self.num_frames = 1;
            self.saving_frequency = 1;
            for step = 1:nSteps
                self.step_size = dt_list(step);
                [~,particles] = self.run_sim("don't save");
                self.initial_particles = particles;
                self.makeImagePsi6Phase(particles, num2str(step) + ".png")
            end
            self.num_frames = old_num_frames;
            self.step_size = old_step_size;
            self.saving_frequency = old_saving_frequency;
            
        end
        
   % ==================================================
   % Utilities/Information
   % ==================================================   
   
        %For use in assessing different starting setups. Returns the
        %fraction of the frame taken up by particles
        %
        %INPUTS:
        %   particles: an nx2 matrix of particles [x y] positions
        %
        %RETURNS:
        %   volFrac: The fraction of the frame taken up by particles
        %
        %
        function volFrac = getVolFrac(self, particles)
            if nargin == 1
                particles = self.initial_particles;
            end
            volFrac = (length(particles)*pi*((self.particle_diam/2)^2))/(self.width*self.height);
        end
        
        
        %Utility for saving pictures and the like. Tells you how many
        %digits the maximum frame number is
        function n = digits_after_t(self)
            n = length(num2str(floor(self.num_frames/self.saving_frequency)));
        end
        
        
        %Tells you the time per frame (in seconds),
        %where frame refers to arrangements saved into plist/cnt_cell 
        %(i.e. we have num_frames/saving_frequency frames)
        function t = time_per_frame(self)
            %Each frame in the simulation is 10 timesteps
            t = 10*self.step_size*self.saving_frequency;
        end
        
        
        
        
   % ==================================================
   % Grain Visualizing
   % ==================================================   
   
   
        %Save a picture of a frame
        %
        %INPUTS
        %   particles: an nx2 matrix of particles [x y] positions
        %   
        %   filename: a string (including .png file extension) for the
        %   image to be saved to.
        %
        %   radius_multiplier (optional): sets how big the particles are dispayed as.
        %   1 is its true size, while something like .6 gives an image that
        %   looks more like our experiments. Defaults to 1 if nothing is
        %   given.
        %
        %OUTPUT
        %   saves an image of the initial setup to the current directory
        %
       function makeImage(self, particles, filename, radius_multiplier)
            
            if nargin == 3
                radius_multiplier = 1;
            end
            
            
            make_dots_wraparound( particles(:,1), particles(:,2), ...
                                  radius_multiplier*self.particle_diam/2, ...
                                  filename, ...
                                  self.width,self.height);
            
       end
        
       %Save a zoomed in picture of a frame. The edges
        %   will be janky, so plan on using a slightly cropped version
        %
        %INPUTS
        %   particles: an nx2 matrix of particles [x y] positions. The
        %   original list. Or a psi6 frame for particles colored by phase.
        %
        %   xlims: [min_x max_x] for your box. Can get this by zooming in
        %   on a figure and typing xlim/ylim on the command line.
        %
        %   ylims: [min_y max_y] for your box
        %
        %   filename: a string (including .png file extension) for the
        %   image to be saved to.
        %   
        %   psi6: make psi6 images or not. Must input psi6_frame instead of
        %   particles in thie case
        %
        %
        %OUTPUT
        %   saves an zoomed in picture to the specified filename. 
        %
       function makeZoomedImage(self, particles, filename, xlims, ylims,psi6)
           if nargin < 6
               psi6 = false;
           end
            
            %get sruff in the box and scale
            in_box = (particles(:,1) > xlims(1) & particles(:,1) < xlims(2) & particles(:,2)...
                < ylims(2) & particles(:,2) > ylims(1));
            particles_in_box = particles(in_box,1:2);
            
            %scale positions
            scaled = zeros(size(particles_in_box));
            scaled(:,1) = (self.width/(xlims(2)-xlims(1))*(particles_in_box(:,1)-xlims(1)));
            scaled(:,2) = (self.height/(ylims(2)-ylims(1))*(particles_in_box(:,2)-ylims(1)));
            
            %find best new radius to draw at
            LC = lattice_constant_fast(scaled, 100);
            radius_multiplier = .98*(LC/self.tightness)/self.particle_diam;
            
            %draw it
            
            if psi6
            not_delta = particles(:,1) ~= 0;
            cmap = generate_thermal_colormap(particles(not_delta,5),-pi,pi, true);
            make_dots_wraparound_color( scaled(:,1), scaled(:,2), ...
                                  radius_multiplier*self.particle_diam/2, ...
                                  filename, ...
                                  self.width,self.height, cmap(in_box,:));
            else
            make_dots_wraparound( scaled(:,1), scaled(:,2), ...
                                  radius_multiplier*self.particle_diam/2, ...
                                  filename, ...
                                  self.width,self.height);    
            end
            
        end
        
        
        %Save a picture of a frame colored by psi6
        %
        %INPUTS
        %   data: either :
        %          psi6_frame: an nx5 matrix [x y psi6mag particle_index psi6phase] 
        %          particle data: an nx2 matrix of particle positions [x y]
        %   savename: a string (including .png file extension) for the
        %   image to be saved to.
        %
        %   radius_multiplier: sets how big the particles are dispayed as.
        %   1 is its true size, while something like .6 gives an image that
        %   looks more like our experiments. Defaults to 1 if nothing is
        %   given.
        %   
        %   use_psi6_sim: a boolean saying whether to use the psi_6_simulation function,
        %                 which gets psi6 for all particles using wraparound boundary conditions,
        %                 instead of cutting off particles at the edge of the image.
        %                 Defaults to false and is only applicable when you give particle_data instead of
        %                 psi6_frame.
        %
        %
        %OUTPUT
        %   saves an image of the initial setup to the current directory
        %
        function makeImagePsi6Phase(self, data, savename,  use_psi6_sim, radius_multiplier)

            %default use_psi6_sim to false
            if nargin < 4
                use_psi6_sim = false;
            end
            
            %check to see if we gave a psi6_frame or frame of a plist
            if size(data,2) == 5
                psi6_frame = data;
            else
                %Make a psi6_frame if it was just position data
                if ~use_psi6_sim
                    TRI = delaunayTriangulation(data(:,1),data(:,2));
                    psi6_frame = psi6_clean_vect(TRI,self.particle_diam*2, self.particle_diam*self.tightness);
                else
                    psi6_frame = psi6_simulation(data(:,1:2),self.width, self.height, self.particle_diam*self.tightness);
                end
            end
            if nargin < 5
                radius_multiplier = 1;
            end
            
            toPlot = psi6_frame(:,1) ~= 0;
            
            cmap = generate_thermal_colormap(psi6_frame(toPlot,5),-pi,pi, true);
            
            make_dots_wraparound_color( psi6_frame(toPlot,1), psi6_frame(toPlot,2), ...
                                  radius_multiplier*self.particle_diam/2, ...
                                  savename, ...
                                  self.width,self.height, cmap);
            
        end
        
        %Save a picture of the initial setup to make sure everything is
        %right
        %
        %INPUTS
        %   savename: a string for the
        %   image to be saved to (excluding .png extension). A normal image
        %   is saved with this name and a copy with p6 phase is saved with
        %   p6 afterwards.
        %
        %   radius_multiplier: sets how big the particles are dispayed as.
        %   1 is its true size, while something like .6 gives an image that
        %   looks more like our experiments. Defaults to 1 if nothing is
        %   given.
        %
        %OUTPUT
        %   saves an image of the initial setup to the current directory
        %
        function preview_setup(self, savename, radius_multiplier)   
            
            %handle optional inputs
            if nargin < 3
                savename = 'initial_setup';
            end
            if nargin < 2
                radius_multiplier = 1;
            end
                
            %make an image
            self.makeImage(self.initial_particles,...
                       strcat(savename,'.png'), radius_multiplier);
            self.makeImagePsi6Phase(self.initial_particles,...
                       strcat(savename,'_psi6.png'),true, radius_multiplier);
        end
        
        
        %Save pictures of a whole experiment
        %
        %INPUTS
        %   cnt_cell: output from the simulations that has been converted
        %   by the pos_to_cnt function.
        %   
        %   savename: a string for the
        %   image to be saved to (no .png). The frame number and .png will be inserted after
        %   the given savename.
        %
        %   radius_multiplier: sets how big the particles are dispayed as.
        %   1 is its true size, while something like .6 gives an image that
        %   looks more like our experiments. Defaults to 1 if nothing is
        %   given.
        %
        %   increment: every x frames in cnt_cell that you want to save
        %   pictures of
        
        %OUTPUT
        %   saves an image of the initial setup to the current directory
        %
        function makeImagesFullExp(self, cnt_cell, savename, increment, radius_multiplier)
            
        if nargin < 4
            increment = 1;
        end
        
        if nargin < 5
            radius_multiplier = 1;
        end

        numFrames = length(cnt_cell);
        
       
        %Go through and get each frame and make a picture for every frame
        for frame = 1:numFrames
            if mod(frame, increment) == 0
                %Get the relevant particles
                inFrame = cnt_cell{frame};
                inFrame = inFrame(:,1:2);

                %Make a filename
                filename = strcat(savename,'_',...
                                      repelem('0',self.digits_after_t()-length(num2str(frame))),...
                                      num2str(frame),'.png');
                
                %Save the image
                self.makeImage(inFrame, filename, radius_multiplier);
                
            end
        end

        end
        
        
        %Save pictures of a whole experiment. Work out of a .mat file for
        %memory savings.
        %
        %INPUTS
        %   cnt_cell_filename: filename of the cnt cell
        %
        %   loading_rate: how many frames to load at a time into memory.
        %   The higher this is the faster it will run, but the greater the
        %   chance you run out of memory. ~500 is a good choice for this
        %   generally.
        %   
        %   savename: a string for the
        %   image to be saved to (no .png). The frame number and .png will be inserted after
        %   the given savename.
        %
        %   radius_multiplier: sets how big the particles are dispayed as.
        %   1 is its true size, while something like .6 gives an image that
        %   looks more like our experiments. Defaults to 1 if nothing is
        %   given.
        %
        %   increment: every x frames in cnt_cell that you want to save
        %   pictures of
        
        %OUTPUT
        %   saves an image of the initial setup to the current directory
        %
        function makeImagesFullExpMAT(self, cnt_cell_filename, loading_rate,...
                                    savename, increment, radius_multiplier)
            
        if nargin < 5
            increment = 1;
        end
        
        if nargin < 6
            radius_multiplier = 1;
        end
        
        %Get the size of the experiment
        m = matfile(cnt_cell_filename);
        [numFrames,~] = size(m,'cnt_cell');
        
        
        %keep track of how much data we have already looked at
        already_processed = 0;

       
        %Go through and get each frame and make a picture for every frame,
        %loading new data when we need it
        for frame = 1:numFrames
            if mod(frame, increment) == 0
                if mod(frame,loading_rate) == 1
                    cnt_cell = m.cnt_cell(frame:min(frame+loading_rate-1,numFrames),1);
                    already_processed = frame - 1;
                end
                %Get the relevant particles
                inFrame = cnt_cell{frame-already_processed};
                inFrame = inFrame(:,1:2);

                %Make a filename
                filename = strcat(savename,'_',...
                                      repelem('0',self.digits_after_t()-length(num2str(frame))),...
                                      num2str(frame),'.png');

                %Save the image
                self.makeImage(inFrame, filename, radius_multiplier);
            end
        end

        end
        
        
        %Save ps6 pictures of a whole experiment
        %
        %INPUTS
        %   psi6_cell: a cell of psi6_frames 
        %       [x y psi6mag particle_index psi6phase]
        %
        %   loading_rate: how many frames to load at a time into memory.
        %   The higher this is the faster it will run, but the greater the
        %   chance you run out of memory. ~500 is a good choice for this
        %   generally.
        %   
        %   savename: a string for the
        %   image to be saved to (no .png). The frame number and .png will be inserted after
        %   the given savename.
        %
        %   radius_multiplier: sets how big the particles are dispayed as.
        %   1 is its true size, while something like .6 gives an image that
        %   looks more like our experiments. Defaults to 1 if nothing is
        %   given.
        %
        %   increment: every x frames in psi6_cell that you want to save
        %   images for
        %OUTPUT
        %   saves an image of the initial setup to the current directory
        %
        function makePsi6ImagesFullExp(self, psi6_cell,...
                                    savename, increment, radius_multiplier)
        
            if nargin < 4
                increment = 1;
            end

            if nargin < 5
                radius_multiplier = 1;
            end

            %find the first index of each frame in pos
            numFrames = length(psi6_cell);


            %Now go through and get each frame and make a picture for every
            %frame
            for frame = 1:numFrames
                
                if mod(frame, increment) == 0

                    %Get the relevant info
                    psi6_frame = psi6_cell{frame};

                    %Make a filename
                    filename = strcat(savename,'_',...
                                      repelem('0',self.digits_after_t()-length(num2str(frame))),...
                                      num2str(frame),'.png');

                    %Save the image
                    self.makeImagePsi6Phase(psi6_frame, ...
                                      filename,false, ...
                                      radius_multiplier);

                end
            end
        end
        
        
        %Save ps6 pictures of a whole experiment
        %
        %INPUTS
        %   psi6_cell: a cell of psi6_frames 
        %       [x y psi6mag particle_index psi6phase]
        %   
        %   savename: a string for the
        %   image to be saved to (no .png). The frame number and .png will be inserted after
        %   the given savename.
        %
        %   radius_multiplier: sets how big the particles are dispayed as.
        %   1 is its true size, while something like .6 gives an image that
        %   looks more like our experiments. Defaults to 1 if nothing is
        %   given.
        %
        %   increment: every x frames in psi6_cell that you want to save
        %   images for
        %OUTPUT
        %   saves an image of the initial setup to the current directory
        %
        function makePsi6ImagesFullExpMAT(self, psi6_cell_filename, loading_rate,...
                savename, increment, radius_multiplier)
        
            if nargin < 5
                increment = 1;
            end

            if nargin < 6
                radius_multiplier = 1;
            end

            
            %Get the size of the experiment
            m = matfile(psi6_cell_filename);
            [numFrames,~] = size(m,'psi6_cell');


            %keep track of how much data we have already looked at
            already_processed = 0;


            %Go through and get each frame and make a picture for every frame,
            %loading new data when we need it
            for frame = 1:numFrames
                if mod(frame, increment) == 0
                    if mod(frame,loading_rate) == 1
                        psi6_cell = m.psi6_cell(frame:min(frame+loading_rate-1,numFrames),1);
                        already_processed = frame - 1;
                    end
                        %Get the relevant info
                        psi6_frame = psi6_cell{frame - already_processed};

                        %Make a filename
                        filename = strcat(savename,'_',...
                                          repelem('0',self.digits_after_t()-length(num2str(frame))),...
                                          num2str(frame),'.png');

                        %Save the image
                        self.makeImagePsi6Phase(psi6_frame, ...
                                          filename, ...
                                          radius_multiplier);
                end
            end
        
        end
        
        
        
   % ==================================================
   % Interfacing with Simulation
   % ==================================================   
        
        %Does a simulation, according to the parameters given
        %
        % Can call this by
        %
        %        sim.run_sim()
        %        sim.run_sim(simulation_number)
        %        sim.run_sim(simulation_number, radius_multiplier, savename)
        %
        %
        %INPUTS
        %   progress_root: a string (with or without .mat file extension) for
        %   the simulation progress to be saved to (in the form of a plist and particle object).
        %
        %   You can also tell it not to save progress by entering "don't
        %   save"
        %
        %   savename: If you want to save images of every frame after completing the simulation.
        %   a string (including .png file extension) for the
        %   images to be saved to. Frame number will be appended to the
        %   end.
        %
        %   radius_multiplier: sets how big the particles are dispayed as.
        %   1 is its true size, while something like .6 gives an image that
        %   looks more like our experiments. Defaults to 1 if nothing is
        %   given.
        %
        %   If no radius_multiplier/savename is given, then the simulation
        %   will run without saving images. 
        %
        %   simulation_number: for when you are running multiple
        %   simulations at once. This indicates which simulation you are
        %   on, used to display progress messages. If no number is given,
        %   defaults to 1.
        %   
        %
        %OUTPUT
        %   particles: the last frame of the experiment:
        %   [x y x_vel y_vel x_accel y_accel particleID]
        %
        %   plist: a history of the particle positions in the experiment
        %   (note that frame is a time stamp):
        %   [x y frame particleID]
        %
        %   Saves images if requested
        %
        function [plist, particles] = run_sim(self, progress_root, ...
                simulation_number, radius_multiplier, savename)
            
            if nargin < 2
                progress_root = "";
            end
            
            if nargin < 3
                simulation_number = 1;
            end
            
            %Get the output files ready
            particles_filename = strcat(progress_root,'particles.mat');
            saving_particles = matfile(particles_filename,'writable',true);
            dontSave = strcmp(progress_root, "don't save");
            if ~dontSave
                saving_particles.particles = self.initial_particles;
            end
            
            %Delete any old stuff
            if ~dontSave
                try
                    delete(particles_filename);
                    delete(plist_filname);
                catch
                    disp('Unable to Delete Old Progress');
                end
            end
            

            
            plist_filename = strcat(progress_root,'plist.mat');        
            saving_plist = matfile(plist_filename,'writable',true);

            num_particles = size(self.initial_particles,1);
            plist = zeros(num_particles*floor(self.num_frames/self.saving_frequency), 4);
            
            %Initialize the first frame
            self.current_frame = 0;
            
            disp('This is from MDSim:');
            disp(size(self.initial_particles));
            
            self.current_particles = self.initial_particles;
            current_row = 1;
            
            %make sure you run chunks so that it's a multiple of saving freq
            n_to_sim_at_a_time = self.saving_frequency*round(self.checkpoint_frequency/self.saving_frequency);
            
            while self.current_frame < self.num_frames
                 
                %run simulation 
                frames_to_simulate =  min(n_to_sim_at_a_time, self.num_frames-self.current_frame);
                [new_particles,new_small_plist]= No_Draw_Simulate_gen(num_particles, self.particle_diam, ...
                    self.neighb_threshold, self.step_size, self.width, self.height, ...
                    self.current_particles, simulation_number,frames_to_simulate, ...
                    self.laser_on, self.laser_off, self.laser_origin, self.saving_frequency);
                %Don't save it directly to allow recovery to work in case
                %of a quit
                small_plist = new_small_plist;
                new_small_plist = 0;
                rows_to_add = size(small_plist,1);
                %Adjust so the plist shows the correct frame number
                small_plist(:,3) = small_plist(:,3) + self.current_frame;
                plist(current_row:current_row+rows_to_add-1,:) = small_plist;
                current_row = current_row+rows_to_add;
                %Save the simulation progress
                self.current_frame = self.current_frame + frames_to_simulate;
                self.current_particles = new_particles;
                particles = self.current_particles;
               %save the particles and plist, unless you say not to
               if ~dontSave
                   disp('hi');
                    save_progress(self.current_frame,saving_plist, plist,...
                        saving_particles, self.current_particles,...
                        current_row, n_to_sim_at_a_time);
               
                    %----------Make image of current frame. Don't if we
                    %actually simulated 0 frames.
                    if ~isempty(particles)
                        self.makeImage(particles(:,1:2), strcat(progress_root, repelem('0',length(num2str(self.num_frames))-length(num2str(self.current_frame))),...
                                          num2str(self.current_frame),'.png'));
                        self.makeImagePsi6Phase(particles(:,1:2),strcat(progress_root,'psi6_', repelem('0',length(num2str(self.num_frames))-length(num2str(self.current_frame))),...
                                          num2str(self.current_frame),'.png'), 1, true);    
                    end
                end
            end
            %save images if desired
            if nargin == 5
                cnt_cell = pos_to_cnt_cell(plist);
                self.makeImagesFullExp(cnt_cell,savename,1,radius_multiplier);
            end
            
        end 
        
        %Resumes a simulation, according to the parameters given
        %
        % Can call this by
        %
        %        sim.resume_sim(particles)
        %        sim.resume_sim(particles, simulation_number)
        %        sim.resume_sim(particles, simulation_number, radius_multiplier, savename)
        %
        %
        %INPUTS
        %   progress_root: a string that is the root of the progress files in the directory.
        %                   i.e. the part before the _plist.mat and
        %                   _particles.mat
        %                   
        %                   you can also enter "don't save" to not save
        %                   progress
        %
        %   framenum: what frame are we resuming from. The simulation will
        %   run for self.num_frames - framenum more frames.
        %
        %   savename: If you want to save images of every frame after completing the simulation.
        %   a string (including .png file extension) for the
        %   images to be saved to. Frame number will be appended to the
        %   end.
        %
        %   radius_multiplier: sets how big the particles are dispayed as.
        %   1 is its true size, while something like .6 gives an image that
        %   looks more like our experiments. Defaults to 1 if nothing is
        %   given.
        %
        %   If no radius_multiplier/savename is given, then the simulation
        %   will run without saving images. 
        %
        %   simulation_number: for when you are running multiple
        %   simulations at once. This indicates which simulation you are
        %   on, used to display progress messages. If no number is given,
        %   defaults to 1.
        %   
        %
        %OUTPUT
        %   particles: the last frame of the experiment:
        %   [x y x_vel y_vel x_accel y_accel particleID
        %
        %   plist: a history of the particle positions in the experiment:
        %   [x y frame particleID]
        %
        %   Saves images if requested
        %
        function [plist, particles] = resume_sim(self,progress_root,...
                                simulation_number, ...
                                radius_multiplier, savename)
            if nargin < 2
                progress_root = "";
                disp('Setting progress_root to be blank');
            end
            if nargin < 3
                simulation_number = 1;
                disp('Setting simulation_number = 1');
            end
           
            disp('In resume_sim now!');
            %load the stuff that's alreay there
            particles_filename = strcat(progress_root,'particles.mat');
            plist_filename = strcat(progress_root,'plist.mat');
            
            disp('These are the filenames right before we try to read in the simulations!');
            disp(particles_filename);
            disp(plist_filename);

            %make sure this is avaliable outside of the try
            loading_particles = 0;
            loading_plist = 0;
            
            disp('trying to read in info from old simulations');
            try
                loading_particles = matfile(particles_filename,'Writable',true);
                loading_plist = matfile(plist_filename,'Writable',true);
                self.current_particles = loading_particles.particles;
                old_plist = loading_plist.plist;
                
            catch
                error("Unable to Load Old Progress in sim " + num2str(simulation_number));
            end
            
            saving_particles = loading_particles;
            saving_plist = loading_plist;
    
            self.current_frame = old_plist(end,3);
            
            %make sure you run chunks so that it's a multiple of saving freq
            n_to_sim_at_a_time = max(self.saving_frequency, self.saving_frequency*round(self.checkpoint_frequency/self.saving_frequency));
            
            num_particles = length(self.initial_particles(:,1));
            
            %Initialize the plist to save stuff to
            current_row = size(old_plist(:,3),1) + 1;
            plist = zeros(num_particles*floor(self.num_frames/self.saving_frequency), 4);
            plist(1:current_row-1,:) = old_plist;
            old_plist = 0;
            
            disp('about to run the simulations! Right before while-loop');
            disp(self.current_frame)
            disp(self.num_frames)
            while self.current_frame < self.num_frames
                disp('This is current_frame: ');
                disp(self.current_frame);
                disp('This is num_frames: ');
                disp(self.num_frames);
                %run simulation
                disp('This is n_to_sim_at_a_time: ');
                disp(n_to_sim_at_a_time);
                frames_to_simulate =  min(n_to_sim_at_a_time, self.num_frames-self.current_frame);
                [new_particles,new_small_plist]= No_Draw_Simulate_gen(num_particles, self.particle_diam, ...
                    self.neighb_threshold, self.step_size, self.width, self.height, ...
                    self.current_particles, simulation_number,frames_to_simulate, ...
                    self.laser_on, self.laser_off, self.laser_origin, self.saving_frequency);
                %Don't save it directly to allow recovery to work in case
                %of a quit
                disp('ran the sim!');
                small_plist = new_small_plist;
                new_small_plist = 0;
                rows_to_add = size(small_plist,1);
                %Adjust so the plist shows the correct frame number
                small_plist(:,3) = small_plist(:,3) + self.current_frame;
                plist(current_row:current_row+rows_to_add-1,:) = small_plist;
                current_row = current_row+rows_to_add;
                %Save the simulation progress
                self.current_frame = self.current_frame + frames_to_simulate;
                disp('This is frames to simulate: ');
                disp(frames_to_simulate);
                self.current_particles = new_particles;
                particles = self.current_particles;
                %save the particles and plist, unless you say not to
               if ~strcmp(progress_root, "don't save")
                   disp("saving progress!");
                   % TODO this line commented out by anna
                   % in order to NOT save to particles and plist
                    %save_progress(self.current_frame,saving_plist, plist,...
                    %saving_particles, self.current_particles,...
                    %current_row, n_to_sim_at_a_time);
               
                    %----------Make image of current frame
                    if ~isempty(particles)
                        disp('Making an image!');
                        self.makeImage(particles(:,1:2), strcat(progress_root, repelem('0',length(num2str(self.num_frames))-length(num2str(self.current_frame))),...
                                          num2str(self.current_frame),'.png'));
                        self.makeImagePsi6Phase(particles(:,1:2),strcat(progress_root,'psi6_', repelem('0',length(num2str(self.num_frames))-length(num2str(self.current_frame))),...
                                          num2str(self.current_frame),'.png'), 1, true);    
                    end
                end
            end
            %save images if desired
            if nargin == 5
                disp('Saving image: ');
                disp(savename);
                cnt_cell = pos_to_cnt_cell(plist);
                self.makeImagesFullExp(cnt_cell,savename,1,radius_multiplier);
            end

        end 
        
        
      
        
    end 
end 



%Saves progress in the form of a plist to the matfile
%saving_plist and self.current_particles to saving_particles.
%
%
function save_progress(current_frame, saving_plist, plist, saving_particles,...
    particles, current_row, n_to_sim_at_a_time)
    %recover if it had quit after saving plist, before writing
    %particles   
    not_right_particles = ~isequal(particles(:,1:2),...
        plist(plist(:,3) == current_frame,1:2));  
    if not_right_particles
        disp("plist correced");
        current_frame = current_frame - n_to_sim_at_a_time;
        plist = plist(plist(:,3) <= current_frame,:);
    end

       function actually_save(first_call, current_frame, saving_plist, plist, particles, saving_particles, current_row)
             %This runs when the function quits, either due to finishing or an error.
             %This second call makes sure that we don't bug out mid saving.
             %Note: this will result in us save over the same thing twice if everything
             %works properly, but it will ensure we save once if stuff goes
             %wrong.
            if first_call
                run_on_quit2 = onCleanup(@()actually_save(false, current_frame, saving_plist, plist, particles, saving_particles, current_row));
            end
            disp('in actually_save!');
            disp(['saving frame ', num2str(current_frame)]);
            saving_plist.plist = plist(1:current_row-1,:);
            disp('plist saved');
            saving_particles.particles = particles;
            disp('particles saved');
       end
   
   %This runs when the function quits, either due to finishing or an error.
   run_on_quit = onCleanup(@()actually_save(true, current_frame, saving_plist, plist, particles, saving_particles, current_row));



end

