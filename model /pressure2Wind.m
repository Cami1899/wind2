%Converts windspeed to pressure based on readings
%Note: Since the observed pressure in flight is reversed (ie pitot tube total
%pressure port goes to negative transducer port), this function inverts the
%windspeeds

function [windSpeed] = pressure2Wind(pressure, rho)
    windSpeed = zeros(size(pressure));
    for i=1:length(pressure)
        if pressure(i) < 0
            windSpeed(i) = sqrt(abs(2*pressure(i)/rho(i)));
        else
            windSpeed(i) = -sqrt(abs(2*pressure(i)/rho(i)));
        end
    end

end