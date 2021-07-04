function done = altsplitmatch()

rad = 5;
tightness = 1.01; % lattice constant = rad*tightness*2
myScale = 2;
xmax = 200*myScale;
ymax = 200*myScale;
R_grain = 50*myScale;
buffer = 5*myScale; % increasing buffer zone seems to help w border madness

% define the boundary of the central enclosed grain
theta = 0:.1:2*pi;
x = xmax/2 + R_grain*cos(theta);
y = ymax/2 + R_grain*sin(theta);
gb_loop = [x',y'];

% get lattice positions associated with the center, left, and right grains
center_lattice = make_hcp(rad*tightness,xmax/2,ymax/2,ceil(xmax/(2*rad*tightness)),0);
left_lattice = make_hcp(rad*tightness,xmax/2+2,ymax/2+10,ceil(xmax/(2*rad*tightness)),10*pi/180);
right_lattice = make_hcp(rad*tightness,xmax/2+3,ymax/2+3,ceil(xmax/(2*rad*tightness)),-15*pi/180);

% each row of DT gives 3 rowIDs that form a delauney triangle
DT = delaunay(center_lattice(:,1),center_lattice(:,2)); 
%triplot(DT,center_lattice(:,1),center_lattice(:,2)); % plot! (if you wanna)

% get positions of particles within GB loop
inLoop = inpolygon(center_lattice(:,1),center_lattice(:,2),gb_loop(:,1), gb_loop(:,2));
banana = center_lattice(inLoop,:); % these are the lattice positions we want. it's called banana bc it's going to split :^)

% get positions of new lattice positions within GB loop
leftSpots = zeros(0,2);
leftMama = zeros(0,2);
for ii = 1:size(left_lattice,1)
    % right half, skip
    if left_lattice(ii,1) > xmax/2+buffer
        continue
    end
    
    distToCenter = sqrt((xmax/2-left_lattice(ii,1))^2+(ymax/2-left_lattice(ii,2))^2);
    % waay inside the gb loop, definitely in leftSpots
    if distToCenter < R_grain - 4*rad*tightness
        leftSpots(end+1,:) = left_lattice(ii,:);
        continue
    end
    
    % waay outside the gb loop, definitely in leftMama
    if distToCenter > R_grain + 4*rad*tightness
        leftMama(end+1,:) = left_lattice(ii,:);
        continue
    end
    
    foundOverlap = 0;
    for jj = 1:size(banana,1)
        dist = sqrt((left_lattice(ii,1)-banana(jj,1))^2 + (left_lattice(ii,2)-banana(jj,2))^2);
        if dist < 2*rad
            leftSpots(end+1,:) = left_lattice(ii,:);
            foundOverlap = 1;
            break
        end
    end
    if ~foundOverlap
        leftMama(end+1,:) = left_lattice(ii,:);
    end
end

rightSpots = zeros(0,2);
rightMama = zeros(0,2);
for ii = 1:size(right_lattice,1)
    % left half, skip
    if right_lattice(ii,1) < xmax/2-buffer
        continue
    end
    
    distToCenter = sqrt((xmax/2-right_lattice(ii,1))^2+(ymax/2-right_lattice(ii,2))^2);
    % waay inside the gb loop, definitely in rightSpots
    if distToCenter < R_grain - 4*rad*tightness
        rightSpots(end+1,:) = right_lattice(ii,:);
        continue
    end
    
    % waay outside the gb loop, definitely in rightMama
    if distToCenter > R_grain + 4*rad*tightness
        rightMama(end+1,:) = right_lattice(ii,:);
        continue
    end
    
    foundOverlap = 0;
    for jj = 1:size(banana,1)
        dist = sqrt((right_lattice(ii,1)-banana(jj,1))^2 + (right_lattice(ii,2)-banana(jj,2))^2);
        if dist < 2*rad
            rightSpots(end+1,:) = right_lattice(ii,:);
            foundOverlap = 1;
            break
        end
    end
    if ~foundOverlap
        rightMama(end+1,:) = right_lattice(ii,:);
    end
end

figure
%mkrSize = 1500;
mkrSize = 300;
scatter(banana(:,1),banana(:,2),mkrSize,'.','MarkerEdgeColor','#CF61BA');
hold on;
axis equal;
xlim([xmax/2-R_grain*1.5 xmax/2+R_grain*1.5]);
ylim([ymax/2-R_grain*1.5 ymax/2+R_grain*1.5]);

%scatter(leftMama(:,1),leftMama(:,2),mkrSize,'.','MarkerEdgeColor','#0072BD');
%scatter(rightMama(:,1),rightMama(:,2),mkrSize,'.','MarkerEdgeColor','#77AC30');
%scatter(leftSpots(:,1),leftSpots(:,2),mkrSize,'.','MarkerEdgeColor','#0072BD');
%scatter(rightSpots(:,1),rightSpots(:,2),mkrSize,'.','MarkerEdgeColor','#77AC30');
%return
% for the new lattice positions on the left and right grids, keep track of which ones are filled
% initially 0 for available. assign to 1 when filled
leftFilled = zeros(size(leftSpots,1),1);
rightFilled = zeros(size(rightSpots,1),1);

% this is called a struct array if you wanna google it
% fields: id?? 
%         pos (position)
%         neighbs (list of nearest neighbor IDs, filled in later)
%         side (this particle goes to the left (-1) or right (+1) lattice. initially 0 for unassigned)
%         wt (explain here lol) left, right
points = struct('pos',num2cell(center_lattice,2),'initPart',num2cell(inLoop,2),'side',0,'wts',zeros(2,1));

% assign nearest neighbors
perms = [1,2,3;2,1,3;3,1,2];
for ii = 1:size(points,1) % loop over center lattice positions
     
    % if there's no particle on this lattice position, skip
    if not(points(ii).initPart)
        continue
    end
    
    %id = parts(ii).id;
    % a list to contain all of this particle's neighbor IDs
    % matlab is very unhappy that i change this list's size but idc
    neighbors = []; 
    for j = 1:size(DT,1) % loop over delauney triangles
        for p = 1:3 % loop over 3 possible positions i might be in
            if DT(j,perms(p,1))==ii
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
    points(ii).neighbs = neighbors;
    
end

% time to start the matching process!
final_positions_left = zeros(0,2);
final_positions_right = zeros(0,2);

% first things first, assign obvious matches
% idk if this makes a difference or not lol
close = 2*rad*tightness/10; % within a/10 means very very close
for ii = 1:size(points,1)
    % if there's no particle on this lattice position, skip
    if not(points(ii).initPart)
        continue
    end
    
    % skip buffer zone particles
    if points(ii).pos(1) < xmax/2+buffer && points(ii).pos(1) > xmax/2-buffer
        continue
    end
    
    myPos = points(ii).pos;
    matchMade = 0;
    
    % search leftSpots for a nearby blue lattice site
    for jj = 1:size(leftSpots,1)
        lattPos = leftSpots(jj,:);
        dist = sqrt( (myPos(1)-lattPos(1))^2 + (myPos(2)-lattPos(2))^2 );
        
        if dist < close
            % NOTE that this is supah inefficient bc assign searches thru
            % leftSpots for the nearest available position, which we
            % already basically did here! but im not gonna fix it unless
            % this program ends up being annoyingly slow lol
            assign(ii,-1);
            matchMade = 1;
            break
        end
    end
    
    if matchMade
        continue
    end
    
    % search rightSpots for a nearby green lattice site
    for jj = 1:size(rightSpots,1)
        lattPos = rightSpots(jj,:);
        dist = sqrt( (myPos(1)-lattPos(1))^2 + (myPos(2)-lattPos(2))^2 );

        if dist < close
            % NOTE that this is supah inefficient bc assign searches thru
            % leftSpots for the nearest available position, which we
            % already basically did here! but im not gonna fix it unless
            % this program ends up being annoyingly slow lol
            assign(ii,1);
            break
        end
    end
end

updateWeights();


% now a while-true loop where assign based on weights
% can use success from assign to check if theres rly an available spot
% also you need to make the semicircles overlap around the GB, otherwise
% there's not enough lattice spots rly
anna = 0;
while 1==1
    anna = anna+1;
    %disp([anna,points(3).wt(1),points(3).wt(2)]);
    for ii = 1:size(points,1) % loop over center lattice positions
        
        % if there's no particle on this lattice position, skip
        if not(points(ii).initPart)
            continue
        end

        % if this particle has already been moved, skip
        if points(ii).side ~= 0
            continue
        end
        
        % skip buffer zone particles
        if points(ii).pos(1) < xmax/2+buffer && points(ii).pos(1) > xmax/2-buffer
            continue
        end
        
        %disp([anna,ii,points(ii).neighbs]);
        % if a lot of your neighbors have made up their minds, its time to
        % pick a side. "a lot" = roughly a sixth of your neighbors
        if points(ii).wt(1)+points(ii).wt(2) >= 0.15
            % if most of ii's neighbors went left, ii should also go left
            if points(ii).wt(1) > points(ii).wt(2)
                assign(ii,-1);
            % if most of ii's neighbors went right, ii should also go right
            elseif points(ii).wt(1) < points(ii).wt(2)
                assign(ii,1);
            % in the event of a tie, decide based on particle ID (basically arbitrary) 
            else
                if mod(ii,2) == 0
                    assign(ii,-1);
                else
                    assign(ii,1);
                end
            end
        end
    end

    updateWeights();
    
    % need to loop until every entry with initPart=1 also has side=+/-1
    % excluding the buffer zone ofc
    allDone = 1;
    for ii = 1:size(points,1)
        if points(ii).initPart && not(points(ii).pos(1) < xmax/2+buffer && points(ii).pos(1) > xmax/2-buffer)
            if points(ii).side == 0
                allDone = 0;
                break
            end
        end
    end
    if allDone == 1
        break
    end
    
end

% finally, go thru and get all the buffer zone particles
for ii = 1:size(points,1) % loop over center lattice positions
        
    % if there's no particle on this lattice position, skip
    if not(points(ii).initPart)
        continue
    end

    % if this particle has already been moved, skip
    if points(ii).side ~= 0
        continue
    end
    
    assign(ii,0);
end

% TODO figure out how to put these BELOW all the quivering plots
scatter(final_positions_left(:,1),final_positions_left(:,2),mkrSize,'.','MarkerEdgeColor','#0072BD');
scatter(final_positions_right(:,1),final_positions_right(:,2),mkrSize,'.','MarkerEdgeColor','#77AC30');

%%%%%%%%%%%%%%% helper functions %%%%%%%%%%%%%%%%%%%%%%
function success = assign(pID, whichSide)
    success = 0;
    
    p1 = points(pID).pos;
    
    if whichSide == 0 % 0 = idc, assign me to either lattice
        nearestIndex = findnearest(points(ii).pos,[leftSpots;rightSpots],[leftFilled;rightFilled]);
        if nearestIndex <= size(leftSpots,1)
            points(pID).side = 1;
            leftFilled(nearestIndex) = 1;
            delta = leftSpots(nearestIndex,:) - p1;
            quiver(p1(1),p1(2),delta(1),delta(2),0,'k','MaxHeadSize',1e2,'AutoScaleFactor',1);
            final_positions_left(end+1,:) = leftSpots(nearestIndex,:);
            success = 1;
        else
            points(pID).side = -1;
            nearestIndex = nearestIndex - size(leftSpots,1);
            rightFilled(nearestIndex) = 1;
            delta = rightSpots(nearestIndex,:) - p1;
            quiver(p1(1),p1(2),delta(1),delta(2),0,'k','MaxHeadSize',1e2,'AutoScaleFactor',1);
            final_positions_right(end+1,:) = rightSpots(nearestIndex,:);
            success = 1;
        end
        return
    end
    
    points(pID).side = whichSide; % -1 = left, 1 = right
    if whichSide == -1
        nearestIndex = findnearest(points(ii).pos,leftSpots,leftFilled);
        leftFilled(nearestIndex) = 1;
        delta = leftSpots(nearestIndex,:) - p1;
        quiver(p1(1),p1(2),delta(1),delta(2),0,'k','MaxHeadSize',1e2,'AutoScaleFactor',1);
        final_positions_left(end+1,:) = leftSpots(nearestIndex,:);
        success = 1;
    elseif whichSide == 1
        nearestIndex = findnearest(points(ii).pos,rightSpots,rightFilled);
        rightFilled(nearestIndex) = 1;
        delta = rightSpots(nearestIndex,:) - p1;
        quiver(p1(1),p1(2),delta(1),delta(2),0,'k','MaxHeadSize',1e2,'AutoScaleFactor',1);
        final_positions_right(end+1,:) = rightSpots(nearestIndex,:);
        success = 1;
    end
    
end

function success = updateWeights()
    
    for kk = 1:size(points,1) % loop over center lattice positions
        % if there's no particle on this lattice position, skip
        if not(points(kk).initPart)
            continue
        end
        
        wt_L = 0; % fraction of neighbors gone left
        wt_R = 0; % fraction of neighbors gone right
        numNeighbs = size(points(kk).neighbs,2);
        for neig = 1:numNeighbs
            
            neighbID = points(kk).neighbs(neig);
            neighbSide = points(neighbID).side;
            
            if neighbSide == -1 % left
                wt_L = wt_L + 1.0/numNeighbs;
            elseif neighbSide == 1 % right
                wt_R = wt_R + 1.0/numNeighbs;
            end
        end
        points(kk).wt = [wt_L;wt_R];
        %points(kk).wt = sum([points(points(kk).neighbs).side]);
    end
    success = 1;
end

done = 1;
end