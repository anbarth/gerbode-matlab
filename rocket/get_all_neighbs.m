function neighbor_data = get_all_neighbs(TRI, LC, delta, far_away)
%gets nearest neighbors of all particles
%INPUTS
%TRI: delaunay triangulation of particle positions
%delta: outline at which we exclude particles
%LC: lattice constant (pixels)
%far_away: how many lattice constants is too far to be a nearest neighb?
%OUTPUTS
%neighbor_data: row j contains the indices of the nearest neighbors of
% particle j

% Handle delta:
if nargin < 3 % default
    delta= [50, 50, 1820, 1116]; 
elseif size(delta) == 1 % scalar
    delta = [delta, delta, 1920 - 2*delta, 1216 - 2*delta]; 
else
    if isfloat(delta)==0 % draw_box disabled for this function, so just use default
        delta = [50 50 1820 1116];
    else % vector
    end
end


%handle far away and LC
if nargin < 2
    [LC, ~] = lattice_constant_fast(TRI.Points,delta);
end
if nargin < 4
    far_away = 1.8;
end

% Get the number of particles
nParticles = length(TRI.Points(:,1));

% Initialize neighbor_data
neighbor_data = zeros(nParticles,8);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Step 2: get neighbor data
%%%%%%%%%%%%%%%%%%%%%%%%

inDelta = (TRI.Points(:,1) > delta(1) & TRI.Points(:,1) < delta(1) + delta(3)) ...
            & (TRI.Points(:,2) > delta(2) & TRI.Points(:,2) < delta(2) + delta(4));

for index = 1:nParticles
    
    % Define the position of the particle:
    home = TRI.Points(index,1:2);
    
    % If we're in delta:
    if inDelta(index)
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %                        START FILTERING                         %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        % Grab the nearest neighbors:
        nns = nearest_neighbors(TRI,index);
        
        % Get their distances to the original particle:
        neighbor_dists = sqrt((TRI.Points(nns,1) - home(1)).^2 ...
            + (TRI.Points(nns,2) - home(2)).^2);
        
        % Only keep ones within a certain distance.
        %disp(nns)
        %disp(neighbor_dists < LC*far_away);
        neighbors = nns(neighbor_dists < LC*far_away);
        neighbors = neighbors(arrayfun(@(pNum) inDelta(pNum), neighbors));
        neighbor_data(index,1:length(neighbors)) = neighbors;
    end
end
end
