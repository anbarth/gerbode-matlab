% Script for plotting the average psi6 of tracked particles, by grain
% Written by Cora Payne '23, Inq Soncharoen '23, Maya Martinez '20
% INPUT: (required in workspace) psi6_cell_tracked.mat, grain_rgb.mat

numG = 2;   % number of grains or granules

psi6_tracked_avg_grain = cell(30,1);

psi6_tracked_avg_grain_plotting = zeros(30,5);
for f=11:40
    psi6_current_frame = psi6_cell_tracked{f};
    
    psi6_avg_frame = zeros(64,6);
    psi6_avg_frame(:,1:2) = psi6_current_frame(:,1:2);
    psi6_avg_frame(:,3) = psi6_current_frame(:,4);
    
    psi6_tracked_avg_grain_plotting(f-10,5) = (f-10)*30-30;    % re-index
    for g=1:numG % g=grain (can change to granule)
        part_in_grain = psi6_current_frame(psi6_current_frame(:,6)==g,:);
        avg = sum(part_in_grain(:,5))/size(part_in_grain,1);
        err = std(part_in_grain(:,5))/sqrt(size(part_in_grain,1));
        
        psi6_tracked_avg_grain_plotting(f-10,g) = avg;
        psi6_tracked_avg_grain_plotting(f-10,g+numG) = err;
        for p=1:64
            if psi6_current_frame(p,6) == g
                psi6_avg_frame(p,4) = avg;
                psi6_avg_frame(p,5) = err;
                psi6_avg_frame(p,6) = g;
            end
        end
    end
    psi6_tracked_avg_grain{f-10} = psi6_avg_frame;
end


figure('Name','Average psi6 phase by grain');

for g=1:numG
    errorbar(psi6_tracked_avg_grain_plotting(1:22,2*numG+1), ...
        psi6_tracked_avg_grain_plotting(1:22,g), psi6_tracked_avg_grain_plotting(1:22,g+numG), ...
        'Color', [grain_rgb(g,1) grain_rgb(g,2) grain_rgb(g,3)]);
    hold on;
end

legend('Grain 1', 'Grain 2')
xlabel('Time Elapsed (sec)');
ylabel('Average \Psi_6 phase');


