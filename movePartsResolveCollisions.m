function done = movePartsResolveCollisions(simName,numSteps)

%folderRoot = "r100/";

srcFold = strcat("../python/crystals/",simName,"/");
%nnFold = strcat("../python/crysNeighbs/",simName,"/");

fstart = strcat(srcFold,"0000.csv");
fend = strcat(srcFold,"end.csv");

% header: row 1 is bead radius
% row 2 is vertices of the polygon that bounds the window
header = readmatrix(fstart,'Range','1:2');

initData = readmatrix(fstart);
finalData = readmatrix(fend);

parts_i = initData(:,2:3);
parts_f = finalData(:,2:3);

inGrain = initData(:,4);

totalDisplacements = parts_f - parts_i;
stepSize = totalDisplacements/numSteps;

for j = 1:numSteps-1
    % set filename
    frameName = compose('%04i',j);
    fname = strcat(srcFold,frameName,".csv");
    
    %nnName = strcat(nnFold,compose('%04i',j),"_neighbs.csv");
    
    % write header
    dlmwrite(fname,[header(1,1)]);
    dlmwrite(fname,header(2,:),'-append');
    
    % displace particles
    parts = parts_i + j*stepSize;
    
    % resolve collisions
    disp(strcat("time to resolve collisions for ",fname,"..."));
    parts = resolveCollisions(parts);
    
    myIDs = ( 1:size(parts,1) )';
    dlmwrite(fname,[myIDs,parts,inGrain],'-append');
    
    writeNeighbs(strcat(simName,"/",frameName));
end

done = 1;
end