%% FindPeaksFunction
%   Takes the concentration data and cell distance vector
%   Returns the number of peaks and the index of each peak.
function [NumPks,PksInd] = findPeaks(data, distance,threshold)
    % Check number of inputs.
    % Fill in unset optional values.
    switch nargin
        case 2
            threshold = max(data)*2.5e-4;
    end
    % Look at the difference between the two points
    diff = data - [0 data(1:length(data)-1)];
    
    % If the difference is significant (> threshold), assign +1 for positive 
    % differences and -1 for negative differences.
    % Otherwise, take on the previous trend value.
    if(diff(1)>0)
        diff(1)=1;
    else
        diff(1)=-1;
    end
    
    for i = 2:length(diff)
        if(diff(i) > threshold)
            diff(i)=1;
        elseif(diff(i)<-threshold)
            diff(i)=-1;
        else
            diff(i)=diff(i-1);
        end
    end
    
    % Identify points where the trend switches from +1 to -1.
    pks = zeros(1,length(distance));
    for i = 1:length(distance)-1
        if(diff(i) > 0 && diff(i+1)<0)
            pks(i)=1;
        end
    end

    % Return number of peaks and their indices.
    PksInd = find(pks);
    NumPks = length(PksInd);
end