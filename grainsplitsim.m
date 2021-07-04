function grainsplitsim(simName,R_grain,miso)

% set up folders named simName in grainsplits, crystals, and crysNeighbs
mkdir(strcat("../grainsplits/",simName))
mkdir(strcat("../python/crystals/",simName))
mkdir(strcat("../python/crysNeighbs/",simName))
mkdir(strcat("../grainsplits/",simName,"/snowflakes"))

% input: grain radius, misorientation angle
[init_prethermal,final_prethermal] = bigsplitmatch(simName,R_grain,miso);
% output: unthermalized particle positions for before and after the split
% stored in init_prethermal.csv and final_prethermal.csv
% also, a specs.txt file that includes details abt this simulation

% thermalize both the initial and final particle positions
thermalize(init_prethermal,simName,"0000",50);
thermalize(final_prethermal,simName,"end",50);
% now there should be crystal files for thermalized particle positions in python/crystals/simName
% as well as corresponding neighbor files in python/crysNeighbs

% fill in intermediate particle positions, resolving collisions as needed
movePartsResolveCollisions(simName,10);
% now there should be crystal files for each timestep in python/crystals/simName
% as well as corresponding neighbor files in python/crysNeighbs

disp("all done!");

end