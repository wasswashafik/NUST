function [channel_effect] = getRandomSubcChannel(association_matrix, BSType)

%change exprnd
subc = getSubCarriers(BSType(1));
channel_effect= zeros(length(association_matrix(1,:)),length(association_matrix(:, 1)), subc );
for BStations = 1:length(association_matrix(1, :));
    
    connected_users = find(association_matrix(:, BStations) == 1);
    for USERS = 1:length(connected_users)
        
        channel_effect(BStations,connected_users(USERS),:) = fading(subc, BSType(BStations));%exprnd(1, [1 subc]);
        
   
    
    end
    
    
    
end


end