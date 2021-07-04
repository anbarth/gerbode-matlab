function writeNeighbs(froot)
%froot = strcat(simName,"/",frameName);
fin = strcat("../python/crystals/",froot,".csv");
fout = strcat("../python/crysNeighbs/",froot,"_neighbs.csv");
A = readmatrix(fin);
A = A(:,2:3); % chop off 1st col and 4th col, if it exists
DT = delaunay(A(:,1),A(:,2)); % each row of DT gives 3 rowIDs that form a delauney triangle
%triplot(DT,A(:,1),A(:,2)); % plot!

perms = [1,2,3;2,1,3;3,1,2];
for i = 1:size(A,1) % loop over particles
    % a list to contain all of i's neighbors 
    % matlab is very unhappy that i change this list's size but idc
    neighbors = [i]; 
    for j = 1:size(DT,1) % loop over delauney triangles
        for p = 1:3 % loop over 3 possible positions i might be in
            if DT(j,perms(p,1))==i
                % other 2 vertices on the triangle are your neighbors
                neighb1 = DT(j,perms(p,2));
                neighb2 = DT(j,perms(p,3));
                % add neighbors only if they're not already in the list
                if not(any(neighbors(:) == neighb1))
                    neighbors(end+1) = neighb1;
                end
                if not(any(neighbors(:) == neighb2))
                    neighbors(end+1) = neighb2;
                end
            end
        end
    end
    if i == 1
        dlmwrite(fout, neighbors);
    else
        dlmwrite(fout, neighbors, '-append');
    end
end
end      

