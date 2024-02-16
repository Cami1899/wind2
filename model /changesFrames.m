%Milo Eirew - September 2023
%Script for changing reference frames from BCM to GPS for either flight

bcmData = bcmPimu1;
bcmLabel = 'bcm1';

gondolaData = gondolaPimu1;
gondolaLabel = 'gondola1';

%allModXY is the matrix containing newly modified (frame-adjusted) pairs of
%Roll rates & Yaw Rates

INS2Tether = 15; %degrees;

pimuZRoll = 6;

%FLIP the pitch and yaw rates, because of the different mounting
%orientations

bcmFlipped = [bcmData(:, 4), -bcmData(:, 5)];

if(strcmp(bcmLabel,'bcm2'))
    rotMatrix = [1/sqrt(2), -1/sqrt(2); 1/sqrt(2),  1/sqrt(2)];
else
    rotMatrix = [0.2588190, 0.9659258; -0.9659258,  0.2588190];
end

allModXY = zeros(length(bcmData), 2);

for i = 1:length(bcmData)
    bcmX = bcmFlipped(i, 1);
    bcmY = bcmFlipped(i, 2);
    
    newXY = rotMatrix*[bcmX; bcmY];

    allModXY(i, :) = newXY;
    a = i

end

bcmSyncARoll2 = sync12(CPU2GPS(bcmLabel, bcmData(:, 1)), allModXY(:, 1), CPU2GPS(gondolaLabel, gondolaData(:, 1))); %SyncA = synced and framed-aligned
bcmSyncAPitch2 = sync12(CPU2GPS(bcmLabel, bcmData(:, 1)), allModXY(:, 2), CPU2GPS(gondolaLabel, gondolaData(:, 1))); %SyncA = synced and framed-aligned
bcmSyncAYaw2 = sync12(CPU2GPS(bcmLabel, bcmData(:, 1)), -bcmData(:, 6), CPU2GPS(gondolaLabel, gondolaData(:, 1))); %SyncA = synced and framed-aligned

plot(CPU2GPS(gondolaLabel, gondolaData(:, 1)), bcmSyncARoll2)
hold on
plot(CPU2GPS(gondolaLabel, gondolaData(:, 1)), gondolaData(:,4))

clf

diffRoll2 = bcmSyncARoll2 - gondolaData(:,4);
diffPitch2 = bcmSyncAPitch2 - gondolaData(:,5);
diffYaw2 = bcmSyncAYaw2 - gondolaData(:,6);

plot(CPU2GPS(gondolaLabel, gondolaData(:, 1)), diffYaw2)
absMagRate1 = (diffRates1(:, 4).^2 + diffRates1(:, 3).^2 + diffRates1(:, 2).^2 ).^0.5;
absMagRate2 = (diffRates2(:, 4).^2 + diffRates2(:, 3).^2 + diffRates2(:, 2).^2 ).^0.5;


return
clf
plot(CPU2GPS(bcmLabel, bcmData(:, 1))-CPU2GPS(bcmLabel, bcmData(1, 1)), allModXY(:, 1), '-o', 'MarkerSize', 0.1, 'Color', '#0062ff')
hold all
plot(CPU2GPS(gondolaLabel, gondolaData(:, 1))-CPU2GPS(bcmLabel, bcmData(1, 1)), gondolaData(:, 4), '-o', 'MarkerSize', 0.1, 'Color', '#ff0000')
legend("BCM", "Gondola")
title("Synced Roll Rates vs Time Since Launch, Flight 1 Ascent", 'FontSize', 20)
ylabel("Roll Rate (rad/s)", 'FontSize', 15)
xlabel("Time since launch (s)", 'FontSize', 15)
set(gcf,'Position',[0 0 1500 500])





