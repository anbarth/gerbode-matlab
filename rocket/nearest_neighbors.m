function nns = nearest_neighbors(TRI,index_num)
%   Outputs the nearest neighbors of a particle with the given index number
%       A vector containing the indicies of the given particles neighbors.
%
%   TRI: the Delaunay Triangulation of a list of particles (using
%   delaunayTriangulation)
%   index_num: the index number of the particle in question

% Obtain the list of triangles:
% A 3N x 3 matrix of indicies. Each row represents a different triangle.
tr = TRI.ConnectivityList;

% Cut this down to triangles containing the particle in question:
tris = tr(tr(:,1) == index_num | tr(:,2) == index_num | tr(:,3) == index_num,:);
% Initialize a list to hold the incoming nearest neighbors:
%{
nns = [];
% Fill the list with nearest neighbors:
for partnum = 1:(size(tris,1)*3) 
    newindex = tris(partnum); % Find the index of the particle located in that position
    if newindex ~= index_num && ~ismember(newindex, nns) % if the particle isn't our starting particle
                                    % and isn't already contained within
                                    % the list:
        nns = [nns, newindex]; % put it into the list   
    end 
end
%}
%simpler and faster code that does the Same Thing
nns = tris(:);
nns = nns(nns ~= index_num);
nns = unique(nns);


end

