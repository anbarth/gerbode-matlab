% line 1: particle radius, idk where to get (ill deal w that outside this script)
% line 2: vertices of boundary (ill deal w that outside this script)
% subsequent lines: x,y for each particle in order
% i'll first do 1-64 in hand ID order
% then i'll do the border particles in... some order

% particles in grain
frame26_positions = info_cell{26}(:,1:2);
dlmwrite("peanut26.csv",frame26_positions);
% boundary particles
dlmwrite("peanut26.csv",info_array_26_outer_left(:,1:2),'-append');
dlmwrite("peanut26.csv",info_array_26_outer_right(:,1:2),'-append');

frame27_positions = info_cell{27}(:,1:2);
dlmwrite("peanut27.csv",frame27_positions);
dlmwrite("peanut27.csv",info_array_26_outer_left(:,1:2),'-append');
dlmwrite("peanut27.csv",info_array_26_outer_right(:,1:2),'-append');