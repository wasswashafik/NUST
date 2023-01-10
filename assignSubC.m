function [subc_association_matrix] = assignSubC(channel_effect, association_matrix)
%channel effect is a 3D array of USERS, BASESTATIONS, and SUBCARRIERS,
%incorporating random effects
%This returns a 3D association matrix
subc_association_matrix = 0.*channel_effect;
[NBS, NU, NSUB] = size(channel_effect);
for BS = 1:NBS
   
    
    connected_users = find(association_matrix(:, BS) == 1);
    
    for j = 1:NSUB
    for USERS = 1:length(connected_users)
        
        
       [maximumVal, index] = max(channel_effect(BS, connected_users(USERS), :));
       if (maximumVal == -10)
           break;
       end
        channel_effect(BS, :, index ) = -10;
        subc_association_matrix(BS, connected_users(USERS), index ) = 1;
    end
    end
    
    
end




end