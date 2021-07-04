function parts = resolveCollisions(parts)

radius = 5;
bumpFactor = 0.5;

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
                numColl = numColl+1;
                parts(ii,:) = parts(ii,:)+(0.5-radius*1.006/delta)*differenceVec*bumpFactor;
                parts(jj,:) = parts(jj,:)-(0.5-radius*1.006/delta)*differenceVec*bumpFactor;
                thereAreCollisions = true;
            end
        end
    end
    counter = counter+1;
end
disp(strcat("resolveCollisions took "+num2str(counter)+" cycles"));
end