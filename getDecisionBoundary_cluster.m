function gb_cell_frame = getDecisionBoundary_cluster(locs, labels, numclusters, LC)
%gets gbs by overlaying a grid over an image of labeled particles and
%finding where the particle label switches
%INPUTS
%locs: [x y] particle matrix
%labels: integer labels representing which particle a grain is in (same
% order as locs)
%numclusters: number of grains
%LC: lattice constant in pixels

%OUTPUTS:
%gb_cell_frame: cell containing gbs for that image (each cell entry is one
% gb, each gb is an [x y] matrix of GB pts)


if nargin < 4
    LC = 25*0.63; % default we usually use, but it's not hard to calculate (lattice constant fast)
end

if nargin < 3
    numclusters = 2;
end

num_pts_x = round(1.5*(max(locs(:,1)) - min(locs(:,1)))/LC);
num_pts_y = round(1.5*(max(locs(:,2)) - min(locs(:,2)))/LC);

% Here is the grid range
%uinds, vinds are grid indices of pts
%u,v are locations on image of pts
[uinds, vinds] = meshgrid(1:num_pts_x, 1:num_pts_y);
u = zeros(num_pts_x,1);
v = zeros(num_pts_y,1);
coords = [uinds(:) vinds(:)];
coords_x = (coords(:,1)-1).*((max(locs(:,1)) - min(locs(:,1)))/num_pts_x) + min(locs(:,1));
coords_y = (coords(:,2)-1).*((max(locs(:,2)) - min(locs(:,2)))/num_pts_y) + min(locs(:,2));
coords_pix = [coords_x, coords_y];
z = zeros(num_pts_y, num_pts_x);
k = dsearchn(locs, coords_pix);
% Evaluate z = theta*x over the grid
for pt_ind = 1:length(k)  
    %get location of real data
    loc_ind = k(pt_ind);
    i = coords(pt_ind,1);
    j = coords(pt_ind,2);
    z(j,i) = labels(loc_ind);
    u(i) = coords_pix(pt_ind, 1);
    v(j) = coords_pix(pt_ind, 2);
end

%{
% remove repeated values
[u, inewu, ~] = unique(u);
v = v(utov(inewu));
z = z(:,utov(inewu));
[v, inewv, ~] = unique(v);
u = u(vtou(inewv));
z = z(vtou(inewv), inewv);
%}

%find each gb
numgbs = numclusters - 1;
gb_cell_frame = cell(numgbs,1);
for gb_ind = 1:numgbs
    boundary_val = gb_ind + 0.5;
    gb = contourc(u, v, z, [boundary_val, boundary_val]);    
    gb = gb.';
    %remove first entry (which is contour height, num_vertices)
    gb = gb(2:end, :);
    %remove points that are clearly wrong
    gb = gb(gb(:,1) < 1920 & gb(:,1) > 0 & gb(:,2) > 0 & gb(:,2) < 1216,:);
    gb_cell_frame{gb_ind,1} = gb;
end


end