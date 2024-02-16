numFiles = 26; %Change depending on how many INS files you are reading
time = [];
allTime = [];

pimu = [];
allPimu = [];

pins = [];
allPins = [];

%Data composed of two alternating datasets: PINS2 & PIMU
% Refer to pages 169-172 of INS documentation:
%https://docs.inertialsense.com/user-manual/reference/user_manual_pdf/InertialSenseDocs.pdf
%Do not assume datasets are always alternating, as that may cause an error
%at the intersection of files
%To index the matrix created by this script, add two to the intended index
%from INS documentation (first two columns are CPU time, )
%Last column of data appears to be hardware flags, disregard

%Change depending on the location of the data directory
basename = 'C:\Users\meirew\Documents\MATLAB\dataProcessing\gondolaData1\INS_';
pinsCount = 0;
pimuCount = 0;
for i = 1:numFiles %Iterates through to find the number of data points, then creates a matrix of the correct size to store the data
    filename = strcat(basename, num2str(i-1), '.DAT');
    fid = fopen(filename);
    tline = fgetl(fid);
    count = 1;
    while ischar(tline)
        data = strsplit(tline, {',', '*'});
        if strcmp(char(data(2)), '$PIMU')
            pimuCount = pimuCount + 1;
        else
            pinsCount = pinsCount + 1;
        end  
        tline = fgetl(fid);
    end
    numDatapoints = i
end

allPimu = zeros(pimuCount, 16);
allPins = zeros(pinsCount, 17);

pimuIndex = 1;
pinsIndex = 1;

for i = 1:numFiles
    filename = strcat(basename, num2str(i-1), '.DAT');
    fid = fopen(filename);
    tline = fgetl(fid);
    count = 1;
    while ischar(tline)
        if width(tline) == 3779
            return
        end
        data = strsplit(tline, {',', '*'});
        if strcmp(char(data(2)), '$PIMU')
            allPimu(pimuIndex, :) = str2double(data);
            pimuIndex = pimuIndex + 1;
        else
            allPins(pinsIndex, :) = str2double(data);
            pinsIndex = pinsIndex + 1;
        end
        tline = fgetl(fid);
        count = count + 1;            
    end

    numDatapoints = i
end
%%




plot(allPimu(:, 1)-allPimu(1, 1), allPimu(:, 5), '-o', 'MarkerSize', 0.1, 'Color', '#ff0000')
title("Gondola Y Axis Roll Rate vs Time Since Startup, Flight 2", 'FontSize', 20)
ylabel("Y Axis Roll Rate (rad/s)", 'FontSize', 15)
xlabel("Time Since Startup (s)", 'FontSize', 15)
set(gcf,'Position',[0 0 1500 500])
