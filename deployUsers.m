function [ user_locations, nusers ] = deployUsers( user_intensity, area )

length = sqrt(area);
center = [length/2 , length/2];         % center coordinates of the circle 
nusers = poissrnd(user_intensity*area);
user_locations = length*rand(nusers,2);



end

