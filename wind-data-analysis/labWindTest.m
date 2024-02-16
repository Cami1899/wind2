%Script for analyzing wind from short lab tests

clf
numFiles = 1;
time = [];
sw1 = [];
sw2 = [];
sw3 = [];
lw1 = [];

allsw1 = []
allsw2 = []
allsw3 = []
alllw1 = []
allBaroTime = []


basename = 'C:\Users\meirew\Documents\above1234';

filename = strcat(basename, '.DAT');
fid = fopen(filename);
tline = fgetl(fid);
count = 1;
while ischar(tline)
    data = strsplit(tline, ',');
    
    time = [time str2double(data(2))];
    sw1 = [sw1 str2double(data(4))];
    sw2 = [sw2 str2double(data(8))];
    sw3 = [sw3 str2double(data(12))];
    lw1 = [lw1 str2double(data(16))];

    tline = fgetl(fid);
    count = count + 1;            
end

clf

pressures = [voltage2Pressure(sw1)', voltage2Pressure(sw2)', voltage2Pressure(sw3)', voltage2Pressure(lw1)'];

rhoLab = 1.225*ones(size(sw1));

plot(time-time(1), pressure2Wind(pressures(:, 1), rhoLab))
hold on
plot(time-time(1), pressure2Wind(pressures(:, 2), rhoLab))
plot(time-time(1), pressure2Wind(pressures(:, 3), rhoLab))
plot(time-time(1), pressure2Wind(pressures(:, 4), rhoLab))


legend("Transducer 1", "Transducer 2", "Transducer 3", "Transducer 4")
