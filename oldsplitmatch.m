function done = oldsplitmatch()

rad = 5;
tightness = 1.01; % lattice constant = rad*tightness
myScale = 2;
xmax = 200*myScale;
ymax = 200*myScale;
R_grain = 50*myScale;

% define the boundary of the central enclosed grain
theta = 0:.1:2*pi;
x = xmax/2 + R_grain*cos(theta);
y = ymax/2 + R_grain*sin(theta);
gb_loop = [x',y'];

% define the boundary of the left semicircular half of the enclosed grain
buffer = 5*myScale;
x_midline_L = [xmax/2+buffer, xmax/2+buffer];
y_midline_L = [ymax/2-sqrt(R_grain^2-buffer^2), ymax/2+sqrt(R_grain^2-buffer^2)];
% TODO theta should actually cover a slightly larger range to account for
% the buffer zone, but i think it doesnt matter as long as the buffer is
% small
theta_L = pi/2:.1:3*pi/2; 
x_L = xmax/2 + R_grain*cos(theta_L);
y_L = ymax/2 + R_grain*sin(theta_L);
xLeft = [x_L,x_midline_L];
yLeft = [y_L,y_midline_L];
left_semicircle_boundary = [xLeft',yLeft'];

% define the boundary of the right semicircular half of the enclosed grain
x_midline_R = [xmax/2-buffer, xmax/2-buffer];
y_midline_R = [ymax/2-sqrt(R_grain^2-buffer^2), ymax/2+sqrt(R_grain^2-buffer^2)];
theta_R = pi/2:-.1:-pi/2;
x_R = xmax/2 + R_grain*cos(theta_R);
y_R = ymax/2 + R_grain*sin(theta_R);
xRight = [x_R,x_midline_R];
yRight = [y_R,y_midline_R];
right_semicircle_boundary = [xRight',yRight'];


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
inLeftHalf = inpolygon(left_lattice(:,1),left_lattice(:,2),left_semicircle_boundary(:,1),left_semicircle_boundary(:,2));
leftSpots = left_lattice(inLeftHalf,:);
inRightHalf = inpolygon(right_lattice(:,1),right_lattice(:,2),right_semicircle_boundary(:,1),right_semicircle_boundary(:,2));
rightSpots = right_lattice(inRightHalf,:);

scatter(banana(:,1),banana(:,2),250,'.','MarkerEdgeColor','#CF61BA');
hold on;
scatter(leftSpots(:,1),leftSpots(:,2),250,'.','MarkerEdgeColor','#0072BD');
scatter(rightSpots(:,1),rightSpots(:,2),250,'.','MarkerEdgeColor','#77AC30');

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

% assign edge particles
for ii = 1:size(points,1) % loop over center lattice positions
    % if there's no particle on this lattice position, skip
    if not(points(ii).initPart)
        continue
    end
    % edge particle! assign to left or right
    if any(inLoop(points(ii).neighbs)==0) % has a neighbor outside GB loop
        if points(ii).pos(1) < xmax/2 % left of center
            assign(ii,-1);
        else % right of center
            assign(ii,1);   
        end
    end
end

updateWeights();
return
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
        %disp([anna,ii,points(ii).neighbs]);
        % if a lot of your neighbors have made up their minds, its time to
        % pick a side. "a lot" = roughly a third of your neighbors
        if points(ii).wt(1)+points(ii).wt(2) >= 0.3
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
    allDone = 1;
    for ii = 1:size(points,1)
        if points(ii).initPart
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



%%%%%%%%%%%%%%% helper functions %%%%%%%%%%%%%%%%%%%%%%
function success = assign(pID, whichSide)
    success = 0;
    points(pID).side = whichSide; % -1 = left, 1 = right
    p1 = points(pID).pos;
    if whichSide == -1
        nearestIndex = findnearest(points(ii).pos,leftSpots,leftFilled);
        leftFilled(nearestIndex) = 1;
        delta = leftSpots(nearestIndex,:) - p1;
        quiver(p1(1),p1(2),delta(1),delta(2),0,'k','MaxHeadSize',1e2,'AutoScaleFactor',1);
        success = 1;
    elseif whichSide == 1
        nearestIndex = findnearest(points(ii).pos,rightSpots,rightFilled);
        rightFilled(nearestIndex) = 1;
        delta = rightSpots(nearestIndex,:) - p1;
        quiver(p1(1),p1(2),delta(1),delta(2),0,'k','MaxHeadSize',1e2,'AutoScaleFactor',1);
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