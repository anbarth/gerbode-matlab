im = imread('200222A_ts_t027.png');
d_im = double(im);
b = bpass(d_im, 1, 20);
pk = pkfnd(b,1,20);


x = pk(:,1);
y = pk(:,2);
%scatter(x,y,'b','filled');


cnt = cntrd(b,pk,25);

leftparts_26 = info_array_26_outer_left(:,1:2);
rightparts_26 = info_array_26_outer_right(:,1:2);

boundaryparts_26 = [leftparts_26; rightparts_26];

for ii = 1:size(boundaryparts_26,1)
    myPos = boundaryparts_26(ii,:);
    %myPos = [693.5, 556];
    closestPos = zeros(1,2);
    minDist = Inf;
    for jj = 1:size(cnt,1)
        dist = sqrt((cnt(jj,1)-myPos(1))^2 + (cnt(jj,2)-myPos(2))^2);
        if dist < minDist
            minDist = dist;
            closestPos = cnt(jj,1:2);
        end
    end
    disp([ii+64,closestPos]);
    %dlmwrite("frame27boundaryparts.csv",closestPos, '-append');
end



