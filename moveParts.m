% WARNING!!!
% this script still uses an old crystal file format
% where pID column didnt exist

folderRoot = "../grainsplits/r150/shitty_rotation/";

srcFold = folderRoot;
%srcFold = strcat("../python/crystals/",folderRoot);
%nnFold = strcat("../python/crysNeighbs/",folderRoot);
numSteps = 6; % num steps from start to end

fstart = strcat(srcFold,"0000.csv");
fend = strcat(srcFold,"end.csv");

% header: row 1 is bead radius
% row 2 is vertices of the polygon that bounds the window
header = readmatrix(fstart,'Range','1:2');

initData = readmatrix(fstart);
finalData = readmatrix(fend);

parts_i = initData(:,2:3);
parts_f = finalData(:,2:3);

partIDs = initData(:,1);
inGrain = initData(:,4);

totalDisplacements = parts_f - parts_i;
stepSize = totalDisplacements/numSteps;

left = floor(min(parts(:,1)));
right = ceil(max(parts(:,1)));
bot = floor(min(parts(:,2)));
top = ceil(max(parts(:,2)));

make_dots_wraparound(parts_i(:,1)-left+5*2,parts_i(:,2)-bot+5*2,5,strcat(srcFold,compose('%04i',0),".png"),right-left+5*3,top-bot+5*3);
for j = 1:numSteps-1
    % set filename
    fname = strcat(srcFold,compose('%04i',j),"_prethermal",".csv");
    %nnName = strcat(nnFold,compose('%04i',j),"_neighbs.csv");
    
    % write header
    dlmwrite(fname,[header(1,1)]);
    dlmwrite(fname,header(2,:),'-append');
    
    % displace particles and record positions
    parts = parts_i + j*stepSize;
    dlmwrite(fname,[partIDs,parts,inGrain],'-append');
    make_dots_wraparound(parts(:,1)-left+5*2,parts(:,2)-bot+5*2,5,strcat(srcFold,compose('%04i',j),"_prethermal",".png"),right-left+5*3,top-bot+5*3);
    %writeNeighbs(fname,nnName);
end
make_dots_wraparound(parts_f(:,1)-left+5*2,parts_f(:,2)-bot+5*2,5,strcat(srcFold,"end",".png"),right-left+5*3,top-bot+5*3);