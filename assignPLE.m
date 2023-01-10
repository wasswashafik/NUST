function [ PLE ] = assignPLE( tier , distance)

PLE = zeros([1 length(tier)]);
for i = 1 : length(tier)
    if tier(i) == 0
        PLE(i) = 2;
    else if tier(i) == 1;
            if isLOS(distance)
                PLE(i) = 2;
            else 
                PLE(i) = 3.3;
            end
        else if tier(i) == 2;
                PLE(i) = 2;
            end
        end
    end

end

