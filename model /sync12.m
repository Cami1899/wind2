%General-purpose function: take in a vector 1 and its associated time, and
%a time 2, and represent the datapoints in vector 1 in the sample points
%given in 2. This is used to compute operations between two vectors which
%both have high sampling frequencies, but sample at different times.

function [returnVec] = sync12(time1, vec1, time2)
    
    returnVec = zeros(size(time2));
    counter1 = 1;

    for i = 1:length(time2)
        if(time2(i) <= time1(1))
            syncedVal = vec1(1);
        elseif (time2(i) >= time1(end))
            syncedVal = vec1(end);
        else
            while(time1(counter1) < time2(i))
                counter1 = counter1 + 1;
            end
            
            syncedVal = mean([vec1(counter1), vec1(counter1-1)]);
        end

        returnVec(i) = syncedVal;

    end
end