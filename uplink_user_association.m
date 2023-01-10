% algorithm to associate each user to a base station defined by its tier.

Ux = (U(:, 1));
Uy = (U(:, 2));
BSx = BSLocation(:, 1);
BSy = BSLocation(:, 2);
do = 1;
uplink_covered_users = 0;
uplink_nU_BS_0 = 0;
uplink_nU_BS_1 = 0;
uplink_nU_BS_2 = 0;
association_matrix_uplink = zeros(length(Ux), length(BSx));


d = [];
Pt_dB = [];
Pr_dB = [];
L_dB = [];
 
for i = 1:length(BSx)
    %This loop runs for each BS
    d(i, :) = sqrt( (Ux - BSx(i)).^2 + (Uy - BSy(i)).^2 ); %distance of users from each BS
   
    Pt_dB = getUplinkTransmitPower(BSType(i)); %get transmit power depending on user type
    L_dB(i, :) = getPathLoss(BSType(i), d(i, :)); %get loss depending on the type of BS connected
    Pr_dB(i, :) = Pt_dB(BSType(i)) - L_dB(i, :);  
    Pr_dB(i, :) = Pr_dB(BSType(i)) + getUplinkAntennaGain(BSType(i), 1); 
  
end

max_power_db = [];
uplink_BS_associated = [];
user_per_tier = [0 0 0]; %index of array -> tier.    value of index -> users associated with that tier
bias_array = [getBiasFactorUplink(0).*ones([1 nBS_0])  getBiasFactorUplink(1).*ones([1 nBS_1])   getBiasFactorUplink(2).*ones([1 nBS_2])];
for i = 1:length(Ux)
    [max_power_db(i), uplink_BS_associated(i)] = max(transpose(Pr_dB(:, i)) + bias_array);
    
  association_matrix_uplink(i, uplink_BS_associated(i)) = 1;
        user_per_tier(1+getTier( uplink_BS_associated(i), nBS_0, nBS_1, nBS_2 )) = user_per_tier(1+getTier( uplink_BS_associated(i), nBS_0, nBS_1, nBS_2 ))+1; 
       
  
    
end
uplink_covered_users = sum(user_per_tier);
uplink_nU_BS_0 = user_per_tier(1);
uplink_nU_BS_1 = user_per_tier(2);
uplink_nU_BS_2 = user_per_tier(3);


