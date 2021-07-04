thetas = [5,10,20];

for ii = 1:length(thetas)
    theta = thetas(ii);
    
    for num = 1:30
        simName = strcat("jul4/r75/theta",theta,"/sim",compose('%02i',num))

        grainsplitsim("jul4/r75/theta5_01"
    end
end