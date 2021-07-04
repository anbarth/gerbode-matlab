function [avgmag, stdmag, psi_6_data] = average_psi6_mag(cnt, delta, ids, picture_name, filetype)
% DESCRIPTION: Outputs the average psi6 magnitude for all the particles in
% a grain as well as the standard deviation of that magnitude and the list
% of psi_6_data for the image.
% INPUTS:
%   cnt:
%   delta (OPTIONAL):
%   ids (OPTIONAL): a list of ids from find_grain_indices
%   picture_name (OPTIONAL):
%   filetype (OPTIONAL):

%OUTPUTS:
% psi_6_data: [DT.Points(index,1) DT.Points(index,2) norm(psi_six_sum) index P];

%The following code was used when this function did not take in cnt
% Prep the image:
%old_cnt = img_prep(picture_name, filetype, focus);
%cnt = filter_dim_particles(old_cnt);

if nargin <2 % delta is not present, so use default
    delta= [50, 50, 1820, 1116]; 
    
elseif size(delta) == 1 % delta is a scalar, so make it into a vector
    delta = [delta, delta, 1920 - 2*delta, 1216 - 2*delta]; 
    
else
    if isfloat(delta)==0 % delta is a string
        delta = draw_box(picture_name, filetype); 

    else % delta is a vector, so do nothing
    end
end

% Gather particle positions:
DT = delaunayTriangulation(cnt(:,1),cnt(:,2));
particle_positions=DT.Points;

% Initialize a list to hold coordinates and psi_6 values:
num_particles=length(particle_positions(:,1));
psi_6_data=zeros(num_particles,5);
psi_six_sum=0;


% Initialize a list containing the index of each particle's nearest
% neighbors:
neighbors_list=zeros(num_particles,10);

% Psi6 code repurposed from Jeremy Wang '17 and also Paul Jerger:
for index = 1:length(DT.Points(:,1))
    if DT.Points(index,1) > delta(1) && DT.Points(index,1) < delta(1) + delta(3)
        if DT.Points(index,2) > delta(2) && DT.Points(index,2) < delta(2) + delta(4)
            
            t = vertexAttachments(DT, index); % triangles including particle "ind"
            t = t{1}; %gets rid of cell format for a list that contains the index of neighboring particles
            N=length(t);
            neighbors = zeros(N,1);
            
            for k = 1:N
                tri = DT.ConnectivityList(t(1,k),:);
                for vert = 1:3 % three vertices to a triangle
                    if tri(1,vert) ~= index
                        for m = 1:N
                            if neighbors(m) == tri(1,vert)
                                break
                            elseif neighbors(m) == 0
                                neighbors(m) = tri(1,vert);
                                neighbors_list(index, m)= tri(1,vert);
                                break
                            end
                        end
                    end
                end
            end
 
            for j = 1:N
                j_vector = DT.Points(neighbors(j),:)- DT.Points(index,:);
                theta =atan2(j_vector(2),j_vector(1));
                psi_six_sum = psi_six_sum + exp(6*theta*sqrt(-1))/N;
            end
            
            % Sets the psi-6 data:
            P=angle(psi_six_sum);
            psi_6_data(index,:)= [DT.Points(index,1) DT.Points(index,2) norm(psi_six_sum) index P];
            psi_six_sum = 0;
        end 
    end
end
imaginary_parts = psi_6_data(psi_6_data(:,1) == 0 & psi_6_data(:,2) == 0 &...
    psi_6_data(:,3) == 0 & psi_6_data(:,4) == 0 & psi_6_data(:,5) == 0,:);
actual_p6d = setdiff(psi_6_data,imaginary_parts,'rows');

if nargin < 3 % If no id list is input, use all of the ids
    avgmag = mean(actual_p6d(:,3));
    stdmag = std(actual_p6d(:,3));
elseif isfloat(ids) == 0
    avgmag = mean(actual_p6d(:,3));
    stdmag = std(actual_p6d(:,3));
elseif size(ids) == 1
    avgmag = mean(actual_p6d(:,3));
    stdmag = std(actual_p6d(:,3));
else
    % Get list of particle positions in the grain:
    grainpos = zeros(length(ids),2);
    for idnum = 1:length(ids)
        id = ids(idnum);
        partpos = cnt(id,1:2);
        grainpos(idnum,:) = partpos;
    end
    
    % Filter actual_p6d to include only those particles:
    new_p6d = actual_p6d(ismember(actual_p6d(:,1:2),grainpos,'rows'),:);
    avgmag = mean(new_p6d(:,3));
    stdmag = std(new_p6d(:,3));
end
           
end