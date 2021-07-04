function done = thermalizeShowAllSteps()
done = 0;
filename = "../grainsplits/r60_theta25/0000_prethermal.csv";
folder = "../grainsplits/r60_theta25/therm_init_smallersteps/";

header = readmatrix(filename,'Range','1:2');
particleData = readmatrix(filename);

parts = particleData(:,2:3); % positions
grainID = particleData(:,4); % 1 for central grain, 0 for mama grains

radius = 5;
makeDotsRad = 5; % make_dots acts like a baby if i give it non-integers

step_size = 0.0025; % typical MDSim value is 0.005
viscosity=1.99*10^-3;
boltz_const = 1.381*(10^(-23));
temp=300;
scaling_factor=2*radius/(1.2*(10^-6)); % px/meter
squig=6*pi*viscosity*(radius/scaling_factor); % 6 pi eta r

% tightness: lattice constant = rad*tightness*2
% only relevant for deciding the size of buffer below
% so its ok if this value is slightly different from in bigsplitmatch.m
tightness = 1.05;
left = floor(min(parts(:,1)));
right = ceil(max(parts(:,1)));
bot = floor(min(parts(:,2)));
top = ceil(max(parts(:,2)));

make_dots_wraparound(parts(:,1)-left+makeDotsRad*2,parts(:,2)-bot+makeDotsRad*2,makeDotsRad,strcat(folder,"0000",".png"),right-left+makeDotsRad*3,top-bot+makeDotsRad*3);
%make_dots_wraparound(parts(:,1),parts(:,2),radius,strcat("thermalize/","0000",".png"),right*2,top*2);

isEdgeParticle = false(size(parts,1),1);
buffer = radius*tightness*10;
for ii = 1:size(parts,1)
    pos = parts(ii,:);
    if pos(1) > right-buffer || pos(1) < left+buffer
        isEdgeParticle(ii) = true;
    elseif pos(2) > top-buffer || pos(2) < bot+buffer
        isEdgeParticle(ii) = true;
    end
end

bumpFactor = 0.5;
for stepNum = 1:25
    random_displacements = normrnd(0,sqrt(2*(boltz_const*temp)/squig*step_size),[size(parts,1) 2])*scaling_factor;
    random_displacements(isEdgeParticle,:)=0;
    parts = parts + random_displacements;

    thereAreCollisions = true;
    counter=1;
    while thereAreCollisions
        thereAreCollisions = false;
        numColl = 0;
        for ii = 1:size(parts,1)
            for jj = ii+1:size(parts,1)

                differenceVec = parts(jj,:)-parts(ii,:);
                delta = sqrt(differenceVec(1)^2+differenceVec(2)^2);

                % collision! resolve!
                if delta < 2*radius*1.005
                    % move each particle back equally until not overlapping
                    % overdo it slightly -- the 1.01 helps avoid floating
                    % point errors
                    numColl = numColl+1;
                    if ~isEdgeParticle(ii)
                        parts(ii,:) = parts(ii,:)+(0.5-radius*1.006/delta)*differenceVec*bumpFactor;
                    end
                    if ~isEdgeParticle(jj)
                        parts(jj,:) = parts(jj,:)-(0.5-radius*1.006/delta)*differenceVec*bumpFactor;
                    end
                    if isEdgeParticle(ii) && isEdgeParticle(jj)
                        parts(ii,:) = parts(ii,:)+(0.5-radius*1.006/delta)*differenceVec*bumpFactor;
                        parts(jj,:) = parts(jj,:)-(0.5-radius*1.006/delta)*differenceVec*bumpFactor;
                    end

                    thereAreCollisions = true;
                end
            end
        end
        %disp(numColl);
        %if counter == 100
        %    break
        %end
        counter = counter+1;
    end
    dlmwrite(strcat(folder,compose('%04i',stepNum),".csv"),header);
    myIDs = ( 1:size(parts,1) )';
    dlmwrite(strcat(folder,compose('%04i',stepNum),".csv"),[myIDs,parts,grainID],'-append');
    make_dots_wraparound(parts(:,1)-left+makeDotsRad*2,parts(:,2)-bot+makeDotsRad*2,makeDotsRad,strcat(folder,compose('%04i',stepNum),".png"),right-left+makeDotsRad*3,top-bot+makeDotsRad*3);
end
done = 1;


end