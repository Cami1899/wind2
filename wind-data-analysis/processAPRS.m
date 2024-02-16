%Process APRS Data

sundayMorning = datetime("2023-07-23 00:00:00");
aprsTime = zeros(1, 94);
aprsAlt = zeros(1, 94);

aprsTemp = zeros(1, 94);
aprsPressure = zeros(1, 94);




for i = 1:height(aprsData)
    timeDiff = seconds(aprsData{i, 1} - sundayMorning);
    aprsTime(i) = timeDiff;

    aprsAlt(i) = aprsData{i, 7};
    message = aprsData{i, 8};
    splitMessage = strsplit(message, {' ', 'C', 'P'});
   

    aprsTemp(i) = str2double(splitMessage(6));
    aprsPressure(i) = str2double(splitMessage(7));
end

%{
plot(aprsTime, aprsAlt, '-o', 'MarkerSize', 1, 'Color', '#ff0000')
title("APRS Altitude vs GPS Time, Flight 2")
xlabel("GPS Time (s)")
ylabel("Altitude (m)")
set(gcf,'Position',[0 0 1500 500])
%}

plot(aprsTime-aprsTime(1), aprsAlt/1000, '-o', 'MarkerSize', 1, 'Color', '#0062ff')
title("APRS Altitude vs Time Since Startup, Flight 2")
xlabel("Time Since Startup (s)")
ylabel("Altitude (m)")
set(gcf,'Position',[0 0 1500 500])


