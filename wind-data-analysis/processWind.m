%Process Wind Data from Flights

clf
load('aprsAlt2.mat')
load('flight2Wind.mat')
load('gondolaAlt1.mat')
load('pressures1.mat')
load('pressures2.mat')

%% 
clc
windArray = flight2Wind; %Voltages
allTime = CPU2GPS('gondola2', windArray(:, 1));


%Pressures1 & Pressures2 are matrices of time, and dynamic pressures
%measures at each of the 4 pitot tubes: SW1`, SW2, SW3, LW1
diff13Pressure1 = pressures1(:, 2)-pressures1(:, 4);
diff24Pressure1 = pressures1(:, 3)-pressures1(:, 5);

diff13Pressure2 = pressures2(:, 2)-pressures2(:, 4);
diff24Pressure2 = pressures2(:, 3)-pressures2(:, 5);

[T1, A1, P1, Rho1] = atmosisa(gondolaAlt1);
[T2, A2, P2, Rho2] = atmosisa(aprsAlt);

syncedRho1 = sync12(CPU2GPS('gondola1', gondolaPins1(1:end-1, 1)), Rho1, CPU2GPS('gondola1', pressures1(:, 1)));
syncedRho2 = sync12(aprsTime, Rho2, CPU2GPS('gondola2', pressures2(:, 1)));


diff13Wind1 = pressure2Wind(diff13Pressure1, syncedRho1);
    
diff24Wind1 = pressure2Wind(diff24Pressure1, syncedRho1);

diff13Wind2 = pressure2Wind(diff13Pressure2, syncedRho2);
    
diff24Wind2 = pressure2Wind(diff24Pressure2, syncedRho2);



windMagnitude1 = (diff13Wind1.^2 + diff24Wind1.^2).^0.5;
windMagnitude2 = (diff13Wind2.^2 + diff24Wind2.^2).^0.5;

dpMagnitude1 = hypot(diff13Pressure1, diff24Pressure1);
dpMagnitude2 = hypot(diff13Pressure2, diff24Pressure2);


plot(CPU2GPS('gondola2', pressures2(:, 1)), windMagnitude2)
absoluteCombinedRoll = (diff13Wind1.^2 + diff24Wind1.^2).^0.5;

syncWindMag2 = sync12(CPU2GPS('gondola2', pressures2(:, 1)), windMagnitude2, CPU2GPS('gondola2', gondolaPimu2(:, 1)));
syncWindMag1 = sync12(CPU2GPS('gondola1', pressures1(:, 1)), windMagnitude1, CPU2GPS('gondola1', gondolaPimu1(:, 1)));

syncPMag2 = sync12(CPU2GPS('gondola2', pressures2(:, 1)), dpMagnitude2, CPU2GPS('gondola2', gondolaPimu2(:, 1)));
syncPMag1 = sync12(CPU2GPS('gondola1', pressures1(:, 1)), dpMagnitude1, CPU2GPS('gondola1', gondolaPimu1(:, 1)));


%{
for i = 1:length(windArray)
    for j = 2:5
        if windArray(i, j) < 2.5
            newPressure = -133*abs((windArray(i, j)/5.0 - 0.5))*(windArray(i, j)/2.0 - 1.25)^2;
        else
            newPressure = 133*abs((windArray(i, j)/5.0 - 0.5))*(windArray(i, j)/2.0 - 1.25)^2;
        end
        pressures2(i, j) = newPressure;
    end
end
%}

return
pressures2(:, 1) = windArray(:, 1);
plot(pressures2(:, 1), pressures2(:, 2), '.r', 'MarkerSize', 0.1)


scale = 2;

startPlotTime = 3000;

windData = [pressures1(:, 1) diff13Wind1, diff24Wind1];


for i = 1:length(pressures1)
    timeElapsed = windData(i, 1)-windData(1,1);
    if mod(i, 10) == 0 && timeElapsed >= startPlotTime 
        clf
        invisPlot=compass([0 scale 0 -scale], [scale 0 -scale 0]);
        set(invisPlot,'Visible','off');
        hold on;
        compass(windData(i, 2), windData(i, 3));
        title("Time Elapsed Since Start = " + timeElapsed)
        drawnow
    end
end


return
plot(allTime-allTime(1), flight1GustSpeeds(2, :))
title("Gust Speed 3 vs Time Since Startup, Flight 2", 'FontSize', 20)
ylabel("Gust Speed 3 (m/s)", 'FontSize', 15)
xlabel("Time Since Startup (s)", 'FontSize', 15)
set(gcf,'Position',[0 0 1500 500])

%atmosisa
%axis([4000 4050 0 4])




