
%Converts CPU Time to GPS Time for any of the avionics stacks involved in the two Lucerne Valley
%Flights. Takes stack label ("bcm1", "bcm2", "gondola1", "gondola2") and
%CPU time vector as arguments

function [gpsTime] = CPU2GPS(stack, time)
    if(stack == "bcm1")
        bcmGPSTime = 396648600/1000;
        bcmCPUTime = 1650028755.0730;
        bcmCPU2GPS = bcmGPSTime - bcmCPUTime;
        gpsTime = time + bcmCPU2GPS;
        return;
    
    
    elseif (stack == "gondola1")
        gondGPSTime = 396167000/1000;
        gondCPUTime = 1649949548.0919;
        gondCPU2GPS = gondGPSTime - gondCPUTime;
        gpsTime = time + gondCPU2GPS;
        return;
    
    
    elseif (stack == "bcm2")
        bcmGPSTime = 397963400/1000;
        bcmCPUTime = 1650035861.8063;
        bcmCPU2GPS = bcmGPSTime - bcmCPUTime;
        gpsTime = time + bcmCPU2GPS;
        return;
    
    
    elseif (stack == "gondola2")
        gondGPSTime = 396840600/1000;
        gondCPUTime = 1649956772.3731;
    
        gondCPU2GPS = gondGPSTime - gondCPUTime;
        gpsTime = time + gondCPU2GPS;
        return;
    

    elseif (stack == "gondolaDecoupled1")
        gondGPSTime = 135284600/1000;
        gondCPUTime = 1625354271.6357;
    
        gondCPU2GPS = gondGPSTime - gondCPUTime;
        gpsTime = time + gondCPU2GPS;
        return;
    

    elseif (stack == "bcmDecoupled1")
        bcmGPSTime = 135172000/1000;
        bcmCPUTime = 1625548671.4032;
    
        bcmCPU2GPS = bcmGPSTime - bcmCPUTime;
        gpsTime = time + bcmCPU2GPS;
        return;
    end
    
    gpsTime = -1

end