function [initFile,finalFile] = bigsplitmatch(folderRoot,R_grain,miso)

initFile = strcat("../grainsplits/",folderRoot,"/0000_prethermal.csv");
finalFile = strcat("../grainsplits/",folderRoot,"/end_prethermal.csv");
specsFile = strcat("../grainsplits/",folderRoot,"/specs.txt");
specsTxt = fopen(specsFile,'w');

rad = 5;
tightness = 1.05;
%R_grain = 100;
thetaL = miso;
thetaR = -1*thetaL;
leftOffset = [-1,1];
rightOffset = [1,5];

fprintf(specsTxt,['grain radius: ' num2str(R_grain) '\n']);
fprintf(specsTxt,['tightness: ' num2str(tightness) '\n']);
fprintf(specsTxt,['theta left: ' num2str(thetaL) '\n']);
fprintf(specsTxt,['theta right: ' num2str(thetaR) '\n']);
fprintf(specsTxt,['particle radius: ' num2str(rad) '\n']);
fprintf(specsTxt,['left offset: ' num2str(leftOffset(1)) ' ' num2str(leftOffset(2)) '\n']);
fprintf(specsTxt,['right offset: ' num2str(rightOffset(1)) ' ' num2str(rightOffset(2)) '\n']);
fclose(specsTxt);


xmax = 2*R_grain+200;
ymax = 2*R_grain+200;
buffer = rad*tightness; 

dlmwrite(initFile,rad);
dlmwrite(finalFile,rad);

% define the boundary of the central enclosed grain
theta = 0:.1:2*pi;
x = xmax/2 + R_grain*cos(theta);
y = ymax/2 + R_grain*sin(theta);
gb_loop = [x',y'];

% define the window
window = [xmax/2-R_grain-rad*tightness*20,ymax/2-R_grain-rad*tightness*20; ...
          xmax/2+R_grain+rad*tightness*20,ymax/2-R_grain-rad*tightness*20; ...
          xmax/2+R_grain+rad*tightness*20,ymax/2+R_grain+rad*tightness*20; ...
          xmax/2-R_grain-rad*tightness*20,ymax/2+R_grain+rad*tightness*20];
     
smallwindow = [xmax/2-R_grain-rad*tightness*4,ymax/2-R_grain-rad*tightness*4; ...
          xmax/2+R_grain+rad*tightness*4,ymax/2-R_grain-rad*tightness*4; ...
          xmax/2+R_grain+rad*tightness*4,ymax/2+R_grain+rad*tightness*4; ...
          xmax/2-R_grain-rad*tightness*4,ymax/2+R_grain+rad*tightness*4];
      
windowToWrite = [smallwindow(1,1),smallwindow(1,2),smallwindow(2,1),smallwindow(2,2),smallwindow(3,1),smallwindow(3,2),smallwindow(4,1),smallwindow(4,2)];      
dlmwrite(initFile,windowToWrite,'-append');
dlmwrite(finalFile,windowToWrite,'-append');

% get lattice positions associated with the center, left, and right grains
center_lattice = make_hcp(rad*tightness,xmax/2,ymax/2,ceil(xmax/(2*rad*tightness)),0); 
left_lattice = make_hcp(rad*tightness,xmax/2+leftOffset(1),ymax/2+leftOffset(2),ceil(xmax/(2*rad*tightness)),thetaL*pi/180); 
right_lattice = make_hcp(rad*tightness,xmax/2+rightOffset(1),ymax/2+rightOffset(2),ceil(xmax/(2*rad*tightness)),thetaR*pi/180);

% get positions of particles within GB loop
inLoop = inpolygon(center_lattice(:,1),center_lattice(:,2),gb_loop(:,1), gb_loop(:,2));
banana = center_lattice(inLoop,:); % these are the lattice positions we want. it's called banana bc it's going to split :^)
%banana = banana( abs(banana(:,1)-xmax/2) > buffer , :);
left_lattice = left_lattice(inpolygon(left_lattice(:,1),left_lattice(:,2),window(:,1),window(:,2)),:);
right_lattice = right_lattice(inpolygon(right_lattice(:,1),right_lattice(:,2),window(:,1),window(:,2)),:);


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
    if distToCenter < R_grain - 0.5*rad*tightness
        leftSpots(end+1,:) = left_lattice(ii,:);
        continue
    end
    
    % waay outside the gb loop, definitely in leftMama
    if distToCenter > R_grain + 0.5*rad*tightness && left_lattice(ii,1) <= xmax/2 
        leftMama(end+1,:) = left_lattice(ii,:);
        continue
    end
    %break
    foundOverlap = 0;
    for jj = 1:size(banana,1)
        dist = sqrt((left_lattice(ii,1)-banana(jj,1))^2 + (left_lattice(ii,2)-banana(jj,2))^2);
        if dist < 1.5*rad % cut out particles that are overlapping _a lot_
            leftSpots(end+1,:) = left_lattice(ii,:);
            foundOverlap = 1;
            break
        end
    end
    if ~foundOverlap && left_lattice(ii,1) <= xmax/2 
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
    if distToCenter < R_grain - 0.5*rad*tightness
        rightSpots(end+1,:) = right_lattice(ii,:);
        continue
    end
    
    foundOverlap = 0;
    % before proceeding, make sure not overlapping w left mama particles
    for jj = 1:size(leftMama,1)
        dist = sqrt((right_lattice(ii,1)-leftMama(jj,1))^2 + (right_lattice(ii,2)-leftMama(jj,2))^2);
        if dist < 1.5*rad
            foundOverlap = 1;
            break
        end
    end
    if foundOverlap
        continue
    end
    
    % waay outside the gb loop, definitely in rightMama
    if distToCenter > R_grain + 0.5*rad*tightness %&& right_lattice(ii,1) >= xmax/2
        rightMama(end+1,:) = right_lattice(ii,:);
        continue
    end
    %break
    foundOverlap = 0;
    for jj = 1:size(banana,1)
        dist = sqrt((right_lattice(ii,1)-banana(jj,1))^2 + (right_lattice(ii,2)-banana(jj,2))^2);
        if dist < 1.5*rad
            rightSpots(end+1,:) = right_lattice(ii,:);
            foundOverlap = 1;
            break
        end
    end
    
    if ~foundOverlap && right_lattice(ii,1) >= xmax/2
        rightMama(end+1,:) = right_lattice(ii,:);
    end
end

allSpots = [leftSpots;rightSpots];

%mkrSize = 5000; % scale 1
%mkrSize = 2000;
mkrSize = 300;
%mkrSize = 1500;
scatter(banana(:,1),banana(:,2),mkrSize,'.','MarkerEdgeColor','#CF61BA');
hold on;
axis equal;
xlim([xmax/2-R_grain-rad*tightness*4 xmax/2+R_grain+rad*tightness*4]);
ylim([ymax/2-R_grain-rad*tightness*4 ymax/2+R_grain+rad*tightness*4]);


scatter(leftMama(:,1),leftMama(:,2),mkrSize,'.','MarkerEdgeColor','#92D5F0');
scatter(rightMama(:,1),rightMama(:,2),mkrSize,'.','MarkerEdgeColor','#C1DB8C');

%scatter(leftSpots(:,1),leftSpots(:,2),mkrSize,'.','MarkerEdgeColor','#000000');
%scatter(rightSpots(:,1),rightSpots(:,2),mkrSize,'.','MarkerEdgeColor','#000000');
%scatter(leftSpots(:,1),leftSpots(:,2),mkrSize,'.','MarkerEdgeColor','#0072BD');
%scatter(rightSpots(:,1),rightSpots(:,2),mkrSize,'.','MarkerEdgeColor','#77AC30');

%return


pythonPartIDLeftMama = ( 1:size(leftMama,1) )';
pythonPartIDRightMama = ( size(leftMama,1)+1:size(leftMama,1)+size(rightMama,1) )';
dlmwrite(initFile, [pythonPartIDLeftMama,leftMama,zeros(size(leftMama,1),1)], '-append');
dlmwrite(finalFile, [pythonPartIDLeftMama,leftMama,zeros(size(leftMama,1),1)], '-append');
dlmwrite(initFile, [pythonPartIDRightMama,rightMama,zeros(size(rightMama,1),1)], '-append');
dlmwrite(finalFile, [pythonPartIDRightMama,rightMama,zeros(size(rightMama,1),1)], '-append');

distances = zeros(size(allSpots,1),size(banana,1));
for ii = 1:size(banana,1)
    for jj = 1:size(allSpots,1)
        dist = sqrt( (banana(ii,1)-allSpots(jj,1))^2 + (banana(ii,2)-allSpots(jj,2))^2 );
        if dist > 2*rad*tightness
            dist = Inf;
        end
        distances(jj,ii) = dist;
    end
    [~,sortedLatticeIDs] = sort(distances(:,ii));
    closestLatticeIDs = sortedLatticeIDs(1:4);
    
    for jj = 1:size(allSpots,1)
        if ~ismember(jj,closestLatticeIDs)
            distances(jj,ii) = Inf;
        end
            
    end
end

%disp(distances(3,94));
%disp(banana(94,1));
%disp(abs(banana(94,1)-xmax/2)<buffer*2.5);
%return

pIDs = zeros(size(banana,1),1);
for ii = 1:size(banana,1)
    pIDs(ii) = ii;
end

latticeIDs = zeros(size(allSpots,1),1);
for jj = 1:size(allSpots,1)
    latticeIDs(jj) = jj;
end


% loop until every particle is assigned a spot
counter = 0;
pythonPartID = size(leftMama,1)+size(rightMama,1)+1;
leftRightCutoff = size(leftSpots,1); % cutoff for lattice IDs, not indices
myQuiver = zeros(size(distances,2),4);
while size(distances,2) > 0
    
    counter = counter+1;
    
    % first, check if any lattice sites are in desperate need
    needyLatticeSite = 0;
    for jj = 1:size(distances,1)
        
        % in buffer zone? skip
        pos = allSpots(latticeIDs(jj),:);
        if abs(pos(1)-xmax/2) < buffer 
            continue
        end
        if sqrt((pos(1)-xmax/2)^2+(pos(2)-ymax/2)^2) > R_grain-rad*tightness
            continue
        end
        
        % if there's only 1 nearby particle left, snatch her up
        isParticleNearby = distances(jj,:) < Inf;
        if sum(isParticleNearby) == 1
            needyLatticeSite = 1;
            % find the one remaining nearby particle
            [~,partInd] = max(isParticleNearby);
            lattInd = jj;
            break
        end
    end
    
    if ~needyLatticeSite
        % assign weights
        wts = zeros(size(distances,2),1);
        for ii = 1:size(distances,2)
            wts(ii) = weight(distances(:,ii));
        end
        
        [~,partInd] = max(wts); % max weight = pickiest particle
        [~,lattInd] = min(distances(:,partInd));
    end
    
    pID = pIDs(partInd);
    particlePos = banana(pID,:); 
    
    latticeID = latticeIDs(lattInd);
    latticeSite = allSpots(latticeID,:);
    
    % if you're trying to vanquish a random GB particle to some random place.... just don't
    if (distances(lattInd,partInd) == Inf) && (abs(particlePos(1)-xmax/2)<buffer*2.5)
        distances = [distances(:,1:partInd-1),distances(:,partInd+1:end)];
        pIDs = [pIDs(1:partInd-1);pIDs(partInd+1:end)];
        scatter(particlePos(1),particlePos(2),mkrSize,'.','MarkerEdgeColor',"#8400FF");
        dlmwrite(initFile, [pythonPartID,particlePos(1),particlePos(2),1], '-append');
        dlmwrite(finalFile, [pythonPartID,particlePos(1),particlePos(2),1], '-append');
        pythonPartID = pythonPartID+1;
        continue
    end
    
    %if needyLatticeSite
    %    myColor = "#DE0021";
    %elseif latticeID <= leftRightCutoff 
    if latticeID <= leftRightCutoff 
        myColor = '#0072BD'; % left is blue
    else
        myColor = '#77AC30'; % right is green
    end
    
    scatter(latticeSite(:,1),latticeSite(:,2),mkrSize,'.','MarkerEdgeColor',myColor);
    delta = allSpots(latticeID,:) - particlePos;
    myQuiver(counter,:) = [particlePos(1),particlePos(2),delta(1),delta(2)];
    dlmwrite(initFile, [pythonPartID,particlePos(1),particlePos(2),1], '-append');
    dlmwrite(finalFile, [pythonPartID,latticeSite(1),latticeSite(2),1], '-append');
    pythonPartID = pythonPartID+1;
    %quiver(p1(1),p1(2),delta(1),delta(2),0,'k','MaxHeadSize',1e2,'AutoScaleFactor',1);
    %text(particlePos(1),particlePos(2),num2str(counter),'FontSize',8);
    %text(particlePos(1),particlePos(2),num2str(pID),'FontSize',8);
    %text(allSpots(latticeID,1),allSpots(latticeID,2),num2str(latticeID),'FontSize',8);
    
    distances = [distances(1:lattInd-1,1:partInd-1),distances(1:lattInd-1,partInd+1:end); ...
                distances(lattInd+1:end,1:partInd-1),distances(lattInd+1:end,partInd+1:end)];
    pIDs = [pIDs(1:partInd-1);pIDs(partInd+1:end)];
    latticeIDs = [latticeIDs(1:lattInd-1);latticeIDs(lattInd+1:end)];
    
    % after matching, if lattice site in buffer zone, remove conflicting lattice sites
    if abs(latticeSite(1)-xmax/2) < buffer
        kk = size(distances,1);
        while kk > 0
            thisSite = allSpots(latticeIDs(kk),:);
            distToThisSite = sqrt((latticeSite(1)-thisSite(1))^2+(latticeSite(2)-thisSite(2))^2);
            if distToThisSite < 1.5*rad
                distances = [distances(1:kk-1,:);distances(kk+1:end,:)];
                latticeIDs = [latticeIDs(1:kk-1);latticeIDs(kk+1:end)];
            end
            kk = kk-1;
        end
    end
end

quiver(myQuiver(:,1),myQuiver(:,2),myQuiver(:,3),myQuiver(:,4),0,'k','MaxHeadSize',1e2,'AutoScaleFactor',1);

% take remaining empty lattice sites around the edge and put them in the mama grains
for jj = 1:size(latticeIDs,1)
    unusedSpot = allSpots(latticeIDs(jj),:);
    distToCenter = sqrt((unusedSpot(1)-xmax/2)^2+(unusedSpot(2)-ymax/2)^2);
    if distToCenter > R_grain-2*rad*tightness && distToCenter < R_grain+2*rad*tightness && abs(unusedSpot(1)-xmax/2) > buffer*2.5
        dlmwrite(initFile, [pythonPartID,unusedSpot(1),unusedSpot(2),0], '-append');
        dlmwrite(finalFile, [pythonPartID,unusedSpot(1),unusedSpot(2),0], '-append');
        pythonPartID = pythonPartID+1;
    end
end

saveas(gcf,strcat("../grainsplits/",folderRoot,"/matched.png"));

    
function wt = weight(L)
    if length(L) < 2
        wt = Inf;
        return
    end
    [sortedDistances,~] = sort(L);
    closestDist = sortedDistances(1);
    secondClosestDist = sortedDistances(2);
    
    wt = (secondClosestDist-closestDist)/closestDist;
end


end