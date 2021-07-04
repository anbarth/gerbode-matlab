function [ sphere_locations ] = make_hcp( rajius,x,y ,n,angle)
%Makes a hexagonal lattice of particles
%INPUTS:
%   rajius: slightly larger than the particle radius
%   x: starting x location
%   y: starting y location
%   n: number of hexagonal rings to step out from the initial point
%   angle: the angle from the horizontal axis you want your lattice to be
%   at (in radians)
%OUTPUTS:
%   sphere_locations: [x y] matrix of locations of particles
%
%-- Jeremy Wang '17
% Commented by Nina Brown '19 & Eli Weissler '19 in summer 2017 
% (this was confusing 2 do so hopefully it makes slightly more sense)
% If you are confused we reccommend drawing a hexagonal lattice out to the
% third ring and walking through the code.

%Initialize a matrix [x y] of particle locations
sphere_locations=zeros(2+3*(n*(n-1)),2);

%Initalize a rotation matrix to rotate a vector by 60 deg (pi/3) clockwise
rotation_matrix=[cos(-pi/3) -sin(-pi/3);sin(-pi/3) cos(-pi/3)];
%rotation sum: keep track of how to get from the 'first' point on a ring to
%the jth point
rotation_sum=[0 0];
%translation: the direction to move out from the initial [x y] point to get
%to the 'first' point on a hexagonal ring
translation=2*rajius*[-cos(-pi/3 + angle) -sin(-pi/3 + angle)];

%step out 1:n from the center point (which hexagonal ring are we in)
for i=1:n
    %2*rajius*rotation initially moves you from the first point to the second
    %point in this hexagonal ring. As the loop goes on its definition
    %changes to take you from the jth 
    rotation=[cos(angle) sin(angle)]';

    %go around the hexagonal ring
    for j=1:6*(i-1)
        
        %get particle index
        index=1 + 3*(i-2)*(i-1) + j;
        %get new sphere location: initial center point + step out to ring 
        %   (starting point for that ring) + rotation sum * particle
        %   diameter (the vector that defines how to get from the starting
        %   point to the jth particle)
        sphere_locations(index,:)=[x y] + (i-1)*translation + 2*rajius*rotation_sum;
        
        %add the rotation vector to rotation_sum
        rotation_sum=rotation_sum + rotation';
        
        % every i-1 particles, we get to a new "side" of the hexagonal
        % ring. Rotate rotation so that it takes us to the next particle
        % once again.
        if mod(j,i-1)==0
            %rotate the rotation vector by pi/3
            rotation=rotation_matrix*rotation;
        end
        
        %so we rotate the vector that goes from from [x y] to the 
        %   translation vector by rotation_sum 
        %we add the angle rotation vector to rotation sum every time we
        %   add a new particle, unless mod(j,i-1) = 0, then we also rotate
        %   the rotation vector by pi/3
    end
    
%add a particle in the middle so we don't get the single weird vacancy
sphere_locations(2+3*(n*(n-1)), 1) = x;
sphere_locations(2+3*(n*(n-1)), 2) = y;
    
end

