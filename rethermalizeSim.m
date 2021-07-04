oldSimName = "r60/theta20";
newSimName = 

% thermalize both the initial and final particle positions
thermalize(strcat("../grainsplits/",simName,"/0000_prethermal.csv"),simName,"0000_redux",25);
thermalize(final_prethermal,simName,"end",25);
% now there should be crystal files for thermalized particle positions in python/crystals/simName
% as well as corresponding neighbor files in python/crysNeighbs

% fill in intermediate particle positions, resolving collisions as needed
movePartsResolveCollisions(simName,10);
% now there should be crystal files for each timestep in python/crystals/simName
% as well as corresponding neighbor files in python/crysNeighbs

disp("all done!");