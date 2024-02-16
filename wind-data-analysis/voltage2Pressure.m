%Converts voltage to pressure
%Refer to Sensirion pressure sensor documentation:
%https://sensirion.com/media/documents/68DF0025/6167E542/Sensirion_Differential_Pressure_Datasheet_SDP8xx_Analog.pdf

function [pressure] = voltage2Pressure(voltage)
    pressure = zeros(size(voltage));
    
    for i = 1:length(voltage)
        
        if(voltage(i) < 2.5)
            pressure(i) = -133*abs(voltage(i)/5 - 0.5)*(voltage(i)/2-1.25)^2;
        else
            pressure(i) = 133*abs(voltage(i)/5 - 0.5)*(voltage(i)/2-1.25)^2;
        end
        
    end

    return
end