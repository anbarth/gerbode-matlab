% in: p, a point
%     lattice, a list of points
%     filled, a list that indicates the availability of points in
%                lattice (0 for available, 1 for filled/unavailable)
% out: the index of the nearest AVAILABLE point to p in lattice
function nearestIndex = findnearest(p,lattice,filled)

minDist = Inf; 
nearestIndex = -1;
for ii=1:size(lattice,1)
    % skip over unavailable spots
    if filled(ii)
        continue
    end
    % this is actually distance^2 but it doesnt matter
    myDist = (p(1)-lattice(ii,1))^2 + (p(2)-lattice(ii,2))^2;
    if myDist < minDist
        minDist = myDist;
        nearestIndex = ii;
    end
end

end

