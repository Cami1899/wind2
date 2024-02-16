%Script for processing Altitude
%Milo Eirew August 2023, adapted from code by Ashish Goel

numFiles = 25;
time = [];
alt = [];

allAlt = [];
allTime = [];

basename = 'C:\Users\meirew\Documents\MATLAB\dataProcessing\bcmData1\INS_';
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
        if strcmp(char(data(2)), '$PINS2')
            % Split the line by commas
            
            % Extract the last two entries (latitude and longitude) from the data
            
            time = [time str2double(data(1))];
            altstr = data(end-1);
            alt = [alt str2double(altstr)];
            
        end
        tline = fgetl(fid);
        count = count + 1;            
    end

    allAlt = [allAlt alt];
    alt = [];

    allTime = [allTime time];
    time = [];

    a = i
end
%%

plot(allTime-allTime(1), allAlt/1000, '.r', 'MarkerSize', 0.1)
title("BCM X axis roll rate vs time since startup, Flight 1", 'FontSize', 20)
xlabel("Time Since Startup (s)", 'FontSize', 15)
ylabel("Altitude (km)", 'FontSize', 15)

set(gcf,'Position',[0 0 1500 500])