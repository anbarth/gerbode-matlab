folderRoot = "r150/";

srcFold = strcat("../python/crystals/",folderRoot);
nnFold = strcat("../python/crysNeighbs/",folderRoot);
numSteps = 10; % num steps from start to end

fstart = strcat(srcFold,"0000.csv");
fend = strcat(srcFold,"end.csv");

step_size = 0.005; % typical MDSim value is 0.005
viscosity=1.99*10^-3;
boltz_const = 1.381*(10^(-23));
temp=300;
scaling_factor=2*radius/(1.2*(10^-6)); % px/meter
squig=6*pi*viscosity*(radius/scaling_factor); % 6 pi eta r
bumpFactor = 0.5;
% header: row 1 is bead radius
% row 2 is vertices of the polygon that bounds the window
header = readmatrix(fstart,'Range','1:2');
radius = header(1);

initData = readmatrix(fstart);
finalData = readmatrix(fend);

parts_i = initData(:,2:3);
parts_f = finalData(:,2:3);

inGrain = initData(:,4);

totalDisplacements = parts_f - parts_i;
stepSize = totalDisplacements/numSteps;

for j = 1:numSteps-1
    disp('-----------------------')
    disp(j)
    % set filename
    fname = strcat(srcFold,compose('%04i',j),".csv");
    nnName = strcat(nnFold,compose('%04i',j),"_neighbs.csv");
    
    % write header
    dlmwrite(fname,[header(1,1)]);
    dlmwrite(fname,header(2,:),'-append');
    
    % displace particles and record positions
    parts = parts_i + j*stepSize;
    

    thereAreCollisions = true;
    counter=1;
    while thereAreCollisions
        thereAreCollisions = false;
        numColl = 0;
        myColors = 0.99*ones(size(parts,1),3);
        stuck = zeros(size(parts,1),1);
        for ii = 1:size(parts,1)
            for jj = ii+1:size(parts,1)

                differenceVec = parts(jj,:)-parts(ii,:);
                delta = sqrt(differenceVec(1)^2+differenceVec(2)^2);

                % collision!
                if delta < 2*radius
                    myColors(ii,:) = [1,0,0];
                    myColors(jj,:) = [1,0,0];
                    stuck(ii,:) = 1;
                    stuck(jj,:) = 1;
                end
            end
        end
        %make_dots_wraparound_color(parts(:,1)-left+makeDotsRad*2,parts(:,2)-bot+makeDotsRad*2,makeDotsRad,strcat(folder,compose('%04i',counter-1),".png"),right-left+makeDotsRad*3,top-bot+makeDotsRad*3,myColors);
        %make_dots_wraparound_color(parts(:,1)+60,parts(:,2)+60,5,strcat(folder,compose('%04i',counter-1),".png"),550,550,myColors);
        disp(counter)
        if mod(counter,100)==0
            random_displacements = normrnd(0,sqrt(2*(boltz_const*temp)/squig*step_size),[size(parts,1) 2])*scaling_factor;
            random_displacements(not(stuck),:)=0;
            parts = parts + random_displacements;
            counter = counter+1;
            thereAreCollisions = true;
            continue
        end
        
        for ii = 1:size(parts,1)
            for jj = ii+1:size(parts,1)

                differenceVec = parts(jj,:)-parts(ii,:);
                delta = sqrt(differenceVec(1)^2+differenceVec(2)^2);

                % collision! resolve!
                if delta < 2*radius
                    % move each particle back equally until not overlapping
                    % overdo it slightly -- the 1.01 helps avoid floating
                    % point errors
                    numColl = numColl+1;
                    %if ~isEdgeParticle(ii)
                    %    parts(ii,:) = parts(ii,:)+(0.5-radius*1.005/delta)*differenceVec*bumpFactor;
                    %end
                    %if ~isEdgeParticle(jj)
                    %    parts(jj,:) = parts(jj,:)-(0.5-radius*1.005/delta)*differenceVec*bumpFactor;
                    %end
                    %if isEdgeParticle(ii) && isEdgeParticle(jj)
                        parts(ii,:) = parts(ii,:)+(0.5-radius*1.005/delta)*differenceVec*bumpFactor;
                        parts(jj,:) = parts(jj,:)-(0.5-radius*1.005/delta)*differenceVec*bumpFactor;
                    %end

                    thereAreCollisions = true;
                end
            end
        end
        %if counter >= 500
        %    break
        %end
        counter = counter+1;
    end
    disp(counter);
    myIDs = ( 1:size(parts,1) )';
    dlmwrite(fname,[myIDs,parts,inGrain],'-append');
    disp(fname);
    writeNeighbs(fname,nnName);
end