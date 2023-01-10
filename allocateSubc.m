function [subc_allocation] = allocateSubc(subc_user_array, channels_per_user)
subc_allocation = 0.*subc_bs_user_array;
for BS = 1:length(subc_bs_user_array(:, 1, 1))
   
    for CPU = 1:  channels_per_user(BS)
        
        for U = 1:subc_bs_user_array(1, : , 1)
            
            [A, index] = max(subc_bs_user_array(BS, U, :));
            subc_allocation(BS, U, index) = 1;
            subc_bs_user_array(BS, U, index) = -1;
            
        end
        
        
    end
    
    
end
end