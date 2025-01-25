function [chQuality] = loadChQuality(fileStr)
   chQuality = zeros(1,96);
fileID = fopen(fileStr,'r');
for k=1:2
    tline = fgets(fileID);
end
for k=3:98
    k
    tline = fgets(fileID);
    chQuality(k-2) = str2num(tline(5:6));
end
fclose(fileID);
end