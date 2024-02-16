%Parse thermistor data and produce arrays
%Battery-adjacent thermistor data will be written to allTemp1
%External thermistor data (if used on stack) will be written to
%allTemp2

% Code by Milo Eirew, adapted from code by Ashish Goel
% August 2023
numFiles = 10;
time = [];
temp1 = [];
temp2 = [];

allTemp1 = []
allTemp2 = []
allTime = []

%Change for your directory
basename = 'C:\Users\meirew\Documents\MATLAB\dataProcessing\gondolaData2\tempSensor_';
for i = 1:numFiles
    filename = strcat(basename, num2str(i-1), '.DAT');
    fid = fopen(filename);
    tline = fgetl(fid);
    count = 1;
    while ischar(tline)
        data = strsplit(tline, ',');
        
        time = [time str2double(data(2))];
        temp1 = [temp1 str2double(data(4))];
        temp2 = [temp2 str2double(data(8))];

        tline = fgetl(fid);
        count = count + 1;            
    end

    allTemp1 = [allTemp1 temp1];
    temp1 = [];

    allTemp2 = [allTemp2 temp2];
    temp2 = [];

    allTime = [allTime time];
    time = [];

    a = i
end
%%

clf

plot(allTime, allTemp2, '-o', 'MarkerSize', 0.1, 'Color', '#ff0000')
title("External Temp (C vs time since start (s)")
xlabel("Time since startup (s)")
ylabel("External Temp (C)")