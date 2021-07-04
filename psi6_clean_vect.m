function [filtered_psi6_data, neighbor_data] = psi6_clean_vect(TRI, delta, LC, far_away)
% Written by: Nina Brown '19
% Description: Calculates the psi6 data for the entered TRI and delta. 
% INPUTS:
%   TRI: The delaunayTriangulation information for the given frame.
%   LC (optional): The lattice constant in pixels. If this is not given,
%   then we will calculate it in the function. 
%   far_away (optional): The condition to exclude far away particles, in
%   LC. If nothing is given, this will default to 1.8. This is both to help
%   when blasting, and to prevent particles at the edge of the image from
%   having 20+ neighbors.
%   delta (optional): an optional input allowing the user to specify a
%   region to include. Defaults to [50 50 1820 1116].
% OUPUTS:
%   filtered_psi6_data: a list of the form [x y magnitude index phase] that
%   does not include contributions to magnitude or phase from far-away
%   neighbors
%   neighbor_data: a list of the form [index1 index2 index3 index4 index5 index6...]
%   where each index corresponds to the nearest neighbor particles of the
%   particle represented by the row number. A 0 represents no more neigbors
%   for that particle.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Step 0: Handle optional inputs %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Handle delta:

if nargin < 2 % default
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
if nargin < 3
    [LC, ~] = lattice_constant_fast(TRI.Points,delta);
end
if nargin < 4
    far_away = 1.8;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Step 1: Setup %
%%%%%%%%%%%%%%%%%

% Get the number of particles
nParticles = length(TRI.Points(:,1));

% Initialize filtered_psi6_data:
filtered_psi6_data = zeros(nParticles,5);

%Step 2: get neighbors
neighbor_data = get_all_neighbs(TRI, LC, delta, far_away);

%Step 3: actually compute psi6
%set up vars
num_cols = length(neighbor_data(1,:));
psi6mat = zeros(length(TRI.Points(:,1)), num_cols);
%loop through each column of neighbor data, getting the psi6 info for each
%point
Nlist = zeros(nParticles, 1);
for col = 1:length(neighbor_data(1,:))
    %set neighbor_data to valid indices
    include_logmat = neighbor_data(:,col) ~= 0;
    
    %get j vectors for each point's next neighbor
    newcol = TRI.Points(neighbor_data(include_logmat,col),:) - TRI.Points(include_logmat,:);
    
    %compute theta and psi6 contribution
    thetacol = atan2(newcol(:,2), newcol(:,1));
    psi6col = exp(6*thetacol*1i);
    
    %add to psi6mat
    psi6mat(include_logmat,col) = psi6col;
    
    %update list of each particles number of nns
    Nlist = Nlist + include_logmat;
end



%compute psi6 vector, phase, and mag
psi6sum = sum(psi6mat, 2);
add_logmat = psi6sum ~= 0;
psi6sum(add_logmat) = psi6sum(add_logmat)./Nlist(add_logmat);

psi6phase = angle(psi6sum);
psi6mag = arrayfun(@(x) norm(x),psi6sum); 

%add to psi6 data matrix
ind_list = (1:nParticles).';
filtered_psi6_data(add_logmat,:) = [TRI.Points(add_logmat,:), psi6mag(add_logmat), ind_list(add_logmat), psi6phase(add_logmat)];
end
