function [filtered_psi6_data, neighbor_data] = psi6_simulation(cnt,width, height, LC, far_away)
% Written by: Eli Weissler '19
% Description: Calculates the psi6 data for the entered cnt. Finds psi6 for
% particles on the edges using wraparound boundary conditions. Takes three
% successive psi6 measurements, shifting stuff to the center so we can see
% it
%
% CURRENT ISSUES: For some reason, it gets the correct psi6 for particles
% near the boundary, but it does not get the correct neighbors.
%
% INPUTS:
%   TRI: The delaunayTriangulation information for the given frame.
%   width: the width of the frame
%   height: the height of the frame
%   LC (optional): The lattice constant in pixels. If this is not given,
%   then we will calculate it in the function. 
%   far_away (optional): The condition to exclude far away particles, in
%   LC. If nothing is given, this will default to 1.8. This is both to help
%   when blasting, and to prevent particles at the edge of the image from
%   having 20+ neighbors.
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
cnt = cnt(:,1:2);
TRI = delaunayTriangulation(cnt);
%handle far away and LC
if nargin < 4
    [LC, ~] = lattice_constant_fast(TRI.Points,1);
end
if nargin < 5
    far_away = 1.8;
end

%Calculate it for a couple of different offsets and then combine

[filtered_psi6_data_1, neighbor_data_1] = psi6_clean_vect(TRI, 3*LC, LC, far_away);

shifted = cnt;
shifted(:,1) = mod(shifted(:,1) + width/3, width);
shifted(:,2) = mod(shifted(:,2) + height/3, height);
TRI = delaunayTriangulation(shifted);
[filtered_psi6_data_2, neighbor_data_2] = psi6_clean_vect(TRI, 3*LC, LC, far_away);
filtered_psi6_data_2(:,1) = mod(filtered_psi6_data_2(:,1) - width/3, width);
filtered_psi6_data_2(:,2) = mod(filtered_psi6_data_2(:,2) - height/3, height);

shifted = cnt;
shifted(:,1) = mod(shifted(:,1) - width/3,width);
shifted(:,2) = mod(shifted(:,2) - height/3, height);
TRI = delaunayTriangulation(shifted);
[filtered_psi6_data_3, neighbor_data_3] = psi6_clean_vect(TRI, 3*LC, LC, far_away);
filtered_psi6_data_3(:,1) = mod(filtered_psi6_data_3(:,1) + width/3, width);
filtered_psi6_data_3(:,2) = mod(filtered_psi6_data_3(:,2) + height/3, height);


%Now combine the frames
[filtered_psi6_data, neighbor_data] = combineArrays(filtered_psi6_data_1,...
                     neighbor_data_1, filtered_psi6_data_2, neighbor_data_2);
[filtered_psi6_data, neighbor_data] = combineArrays(filtered_psi6_data,...
                     neighbor_data, filtered_psi6_data_3, neighbor_data_3); 
                 
%make sure stuff is symmetric
for p = 1:size(neighbor_data,1)
    %get all nns
    nns = neighbor_data(p, neighbor_data(p,:)~=0);
    for n = 1:length(nns)
        %if the other particle doesn't have p as a neighbor, add it
        if ~ismember(p,neighbor_data(nns(n),:))
            last_n = find(neighbor_data(nns(n),:),1,'last');
            neighbor_data(nns(n),last_n+1) = p;
        end
    end

end
end
%helper function that combines the data from two different psi6 or
%neighbor_data. It looks at each row and keeps the entries from the one
%with more neighbors, unless the p6_row is empty, in which case it takes
%the non empty one.
function [p6_out, nb_out] = combineArrays(p61,nb1,p62,nb2)
    %make sure everything is the same length
    num_particles = size(p61,1);
    assert(num_particles == size(p62,1));
    
    %initialize the output
    p6_out = zeros(num_particles, 5);
    nb_out = zeros(num_particles, max([size(nb1,2) size(nb2,2)]));
    
    for i = 1:num_particles
        p61_row = p61(i,:);
        nb1_row = nb1(i,:);
        p62_row = p62(i,:);
        nb2_row = nb2(i,:);
                
        %check if one of the psi6 rows is 0. This indicates it was inside
        %of delta for one of the frames, so pick the other.
        if p61_row(3) == 0
            p6_out(i,:) = p62_row;
            nb_out(i,1:length(nb2_row)) = nb2_row;
        elseif p62_row(3) == 0
            p6_out(i,:) = p61_row;
            nb_out(i,1:length(nb1_row)) = nb1_row;
        else
            %if niether are zero, then pick the one that's closer to 6 nns
            nn1 = length(find(nb1_row));
            nn2 = length(find(nb2_row));
            if nn1 == 6
                p6_out(i,:) = p61_row;
                nb_out(i,1:length(nb1_row)) = nb1_row;
            elseif nn2 == 6
                p6_out(i,:) = p62_row;
                nb_out(i,1:length(nb2_row)) = nb2_row;
            elseif abs(nn1-6) < abs(nn2-6)
                p6_out(i,:) = p61_row;
                nb_out(i,1:length(nb1_row)) = nb1_row;
            else
                p6_out(i,:) = p62_row;
                nb_out(i,1:length(nb2_row)) = nb2_row;
            end
        end
            
    end

end







