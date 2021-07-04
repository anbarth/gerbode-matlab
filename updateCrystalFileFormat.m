
folderName = "../python/crystals/bigboy/";
fileList = dir(folderName);

for i=1:size(fileList,1)
    fileName = fileList(i).name;
    disp(fileName);
    
    if strlength(fileName)<5 || fileName(end-3:end) ~= ".csv"
        continue
    end
    if fileName == "splitmatch_init.csv"
        continue
    end
    
    fullFileName = strcat(folderName,fileName);
    header = readmatrix(fullFileName,'Range','1:2');
    positions = readmatrix(fullFileName);
    
    numParts = size(positions,1);
    partIDs = (1:numParts)';
    
    dlmwrite(fullFileName,header(1,1));
    dlmwrite(fullFileName,header(2,:),'-append');
    dlmwrite(fullFileName,[partIDs,positions],'-append');
end
