function [a, std_dev] = lattice_constant_fast(cnt, delta)
% DESCRIPTION: Finds the average lattice constant for an input cnt.
% INPUTS:
%   cnt: can be vector of the form [x y brightness radius_of_gyration] or
%   [x y] or just anything with positions in the first two columns.
%   delta: (optional) allows the user to specify the area the code
%   should be run on with a vector of the form [xmin ymin width height]. If
%   left blank, defaults to [50 50 1820 1116]; if a scalar is entered,
%   defaults to [scalar scalar 1920-2*scalar 1216-2*scalar]

% Step 0: Handle optional inputs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Handle delta:
if nargin < 2
    delta= [50, 50, 1820, 1116]; % default delta
elseif size(delta) == 1 % delta is a scalar
    delta = [delta, delta, 1920 - 2*delta, 1216 - 2*delta]; 
else
    if isfloat(delta)==0 % delta is a string
        disp('Box selection disabled for lattice_constant_fast!');
    else % delta is a vector, so do nothing
    end
end

% Step 1: Generate a list of orderly particles to use for the calculation:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Find psi6:
[~, ~, psi6d] = average_psi6_mag(cnt);

% Get the Delaunay Triangles for all the particles:
TRI = delaunayTriangulation(cnt(:,1),cnt(:,2));

% Generate index_list:
%index_list = zeros(length(cnt(:,1)));
%for index = 1:length(index_list)
    %index_list(index) = index;
%end
nonzero_cnt = cnt(cnt(:,1) ~= 0 ,:);
index_list = 1:length(nonzero_cnt(:,1));
% Filter the indices based on psi6 (we want only the most orderly particles):
n = 1;
while n <= length(index_list)
    index = index_list(n);
    psi6mag = psi6d(index, 3);
    if psi6mag < 0.9
        index_list = [index_list(1:n-1),index_list(n+1:length(index_list))];
    else
        n = n+1;
    end
end

% Step 2: Find the lattice constant:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initialize a list to hold average distances:
avgdistlist = zeros(length(index_list));

% Find the average distances:
for inum = 1:length(index_list)
    partnum = index_list(inum); % Grab the particle index
    home = TRI.Points(partnum,:); % Find the position of that particle
    if home(1) > delta(1) && home(1) < delta(1) + delta(3) && home(2) > delta(2) &&...
            home(2) < delta(2) + delta(4) % the particle is contained in the box
        neighbors = nearest_neighbors(TRI,partnum); % Get that particle's nearest neighbors
        poslist = zeros(size(neighbors,2),2); % Initialize a list to hold the positions
                                            % of the nns
        for n = 1:size(neighbors,2) % Fill the list:
            neighborpos = TRI.Points(neighbors(n),:); 
            poslist(n,:) = neighborpos;
        end
        % Find the distance between the positions in the list and the position
        % of the original particle:
        dist_list=sqrt((poslist(:,1) - home(1)).^2 + (poslist(:,2)- home(2)).^2); 
        avgdist = mean(dist_list); % Average the list of distances
        avgdistlist(inum) = avgdist;
    else
    end
end

% Cut out the zeros:
finaldistlist = avgdistlist(avgdistlist ~= 0);

a = mean(finaldistlist); % Average the list of distances of all the particles;
                        % let this be the lattice constant
a = a(1); % For some reason it spits out a vector?? Who knows why.
std_dev = std(finaldistlist); % Standard deviation for good measure
std_dev = std_dev(1);
end
    
