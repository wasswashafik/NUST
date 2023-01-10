function [interference, interference_dB] = getDownlinkInterference...
    (BS_associated, BSType, nBS_0, nBS_1, nBS_2 ,downlink_association_matrix, Downlink_Pr_matrix, angle_array, beam_width)


%for each user, sum power being received at that user from all BS excluding
%the tagged BS
nU = length(downlink_association_matrix(:, 1));
interference = zeros([1 nU]); %Size same as number of users, interference for each user obtained
power_watts = 10.^(Downlink_Pr_matrix/10);

%add fading to each link

AntennaAlignedGain =  getDownlinkAntennaGain(1, 1);
AntennaMisAlignedGain =  getDownlinkAntennaGain(1, 0);

AntennaAlignedGain = 10.^(AntennaAlignedGain/10);
AntennaMisAlignedGain = 10.^(AntennaMisAlignedGain/10);
for user = 1:nU % for each user
 
    if BSType(BS_associated(user)) == 0
    interference(user) = sum(power_watts(user, (1:nBS_0))) + sum(power_watts(user, (nBS_1+nBS_0+1):(nBS_0+nBS_1+nBS_2))) - power_watts( user , BS_associated(user));
    
    elseif BSType(BS_associated(user)) == 1
        for B_S = (nBS_0+1):(nBS_0+nBS_1)
            if B_S == BS_associated(user)
                continue;
            end
        presentangle = (360.*rand()) - 180;
        
        if ( angle_array(user,BS_associated(user)) < presentangle + (beam_width/2) && angle_array(user,BS_associated(user))> presentangle - (beam_width/2))         
            interference(user) = interference(user) + AntennaAlignedGain.*power_watts(user, B_S);  
        else
            interference(user) = interference(user) + AntennaMisAlignedGain.*power_watts(user, B_S);
        end
        end
        
    elseif BSType(BS_associated(user)) == 2
        
         interference(user) = sum(power_watts(user, (1:nBS_0))) + sum(power_watts(user, (nBS_1+nBS_0+1):(nBS_0+nBS_1+nBS_2))) - power_watts( user , BS_associated(user));
    end
    
    
end


interference_dB = 10.*log10(interference);


end