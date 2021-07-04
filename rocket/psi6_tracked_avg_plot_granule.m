% Script for plotting the average psi6 of tracked particles, by granule
% Written by Cora Payne '23, Inq Soncharoen '23, Maya Martinez '20
% INPUT: (required in workspace) psi6_cell_tracked.mat, granule_rgb.mat

numG = 10;   % number of grains or granules

psi6_tracked_avg_granule = cell(30,1);

psi6_tracked_avg_granule_plotting = zeros(30,2*numG+1);
for f=11:40
    psi6_current_frame = psi6_cell_tracked{f};
    
    psi6_avg_frame = zeros(64,6);
    psi6_avg_frame(:,1:2) = psi6_current_frame(:,1:2);
    psi6_avg_frame(:,3) = psi6_current_frame(:,4);
    
    psi6_tracked_avg_granule_plotting(f-10,2*numG+1) = (f-10)*30-30;    % re-index
    for g=1:numG % g=grain (can change to granule)
        part_in_g = psi6_current_frame(psi6_current_frame(:,7)==g,:);
        avg = sum(part_in_g(:,5))/size(part_in_g,1);
        err = std(part_in_g(:,5))/sqrt(size(part_in_g,1));
        
        psi6_tracked_avg_granule_plotting(f-10,g) = avg;
        psi6_tracked_avg_granule_plotting(f-10,g+numG) = err;
        for p=1:64
            if psi6_current_frame(p,7) == g
                psi6_avg_frame(p,4) = avg;
                psi6_avg_frame(p,5) = err;
                psi6_avg_frame(p,6) = g;
            end
        end
    end
    psi6_tracked_avg_granule{f-10} = psi6_avg_frame;
    
    
end

for g=1:numG
    disp(g);
    disp(psi6_tracked_avg_granule_plotting(16,g));
    disp(psi6_tracked_avg_granule_plotting(17,g));
    disp('---------------------')
end

figure('Name','Average psi6 phase by granule');

for g=1:numG
    %errorbar(psi6_tracked_avg_granule_plotting(1:22,2*numG+1), ...
    %    psi6_tracked_avg_granule_plotting(1:22,g), psi6_tracked_avg_granule_plotting(1:22,g+numG), ...
    %    'Color', [granule_rgb(g,1) granule_rgb(g,2) granule_rgb(g,3)]);
    % cutting out color bc i dont have granule_rgb
    errorbar(psi6_tracked_avg_granule_plotting(1:22,2*numG+1), ...
        psi6_tracked_avg_granule_plotting(1:22,g), psi6_tracked_avg_granule_plotting(1:22,g+numG));
    hold on;
end

legend('Granule 1', 'Granule 2', 'Granule 3', 'Granule 4', 'Granule 5', 'Granule 6', 'Granule 7', 'Granule 8', 'Granule 9', 'Granule 10');
xlabel('Time Elapsed (sec)');
ylabel('Average \Psi_6 phase');