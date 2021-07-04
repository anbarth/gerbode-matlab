R_grain = 60;
tightness = 1.05;
rad = 5;

xmax = 2*R_grain+200;
ymax = 2*R_grain+200;
buffer = rad*tightness; 

theta = 0:.25:2*pi;
x = xmax/2 + (R_grain+2*buffer)*cos(theta);
y = ymax/2 + (R_grain+2*buffer)*sin(theta);
myCirc = [x',y'];

angleNames = ["5","6p25","7p5","8p75","10","15","20","25"];
simNames = strings(1,size(angleNames,2));
for ii = 1:size(angleNames,2)
    simNames(ii) = strcat("r60/theta",angleNames(ii));
end

for ii = 1:size(simNames,2)
    simName = simNames(ii);
    initFile = "../python/crystals/"+simName+"/0000.csv";
    initData = readmatrix(initFile);
    partIDs = initData(:,1);
    parts = initData(:,2:3);
    inCircle = inpolygon(parts(:,1),parts(:,2),myCirc(:,1), myCirc(:,2));
    dlmwrite(strcat("../grainsplits/",simName,"/partsInCircle.csv"),[partIDs,inCircle]);
end
% window = zeros(1,2*size(x,1));
% for ii = 1:size(x,2)
%     window(2*ii-1) = x(ii);
%     window(2*ii) = y(ii);
% end
% 
% dlmwrite("../grainsplits/r60/circularWindow.csv", window);