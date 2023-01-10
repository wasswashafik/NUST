function [angle_array] = getBsToUserAngle(U, nU, nBS_0, nBS_1, BSLocation)

Ux = (U(:, 1));
Uy = (U(:, 2));
BSx = BSLocation(:, 1);
BSy = BSLocation(:, 2);

angle_array = zeros(nU, length(BSx));

for BS = nBS_0+1:(nBS_0+nBS_1)
    
    for user = 1:nU
        
       angle_array (user, BS) = (180./pi).*angle(Ux(user) - BSx(BS) + 1j.*(Uy(user) - BSy(BS)) );
        
    end
    
end

end