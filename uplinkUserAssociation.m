function [uplink_BS_associated, uplink_nU_BS_0, uplink_nU_BS_1, uplink_nU_BS_2, association_matrix_uplink, Pr_dB, max_power_db ] = uplinkUserAssociation(Pt_dB, U, nU, nBS_0, nBS_1, nBS_2, BSType, BSLocation)

Ux = (U(:, 1));
Uy = (U(:, 2));

BSx = BSLocation(:, 1);
BSy = BSLocation(:, 2);

uplink_nU_BS_0 = 0;
uplink_nU_BS_1 = 0;
uplink_nU_BS_2 = 0;

association_matrix_uplink = zeros(length(Ux), length(BSx));
d = zeros([nBS_0+nBS_1+nBS_2 , nU]);
Pr_dB = zeros([nBS_0+nBS_1+nBS_2 ,nU]);
L_dB = zeros([nBS_0+nBS_1+nBS_2 ,nU]);

for i = 1:length(BSx)
    d(i, :) = sqrt( (Ux - BSx(i)).^2 + (Uy - BSy(i)).^2 ); %distance of users from each BS
    L_dB(i, :) = getPathLoss(BSType(i), d(i, :)); %get loss depending on the type of BS connected
    Pr_dB(i, :) = Pt_dB(1, :) - L_dB(i, :);
    Pr_dB(i, :) = Pr_dB(i, :) + getUplinkAntennaGain(BSType(i), 1);
end

max_power_db = [];
uplink_BS_associated = [];
user_per_tier = [0 0 0]; %index of array -> tier.    value of index -> users associated with that tier
bias_array = [getBiasFactorUplink(0).*ones([1 nBS_0])  getBiasFactorUplink(1).*ones([1 nBS_1])   getBiasFactorUplink(2).*ones([1 nBS_2])];

for i = 1:nU
    [max_power_db(i), uplink_BS_associated(i)] = max(transpose(Pr_dB(:, i)) + bias_array);
    association_matrix_uplink(i, uplink_BS_associated(i)) = 1;
    user_per_tier(1+getTier( uplink_BS_associated(i), nBS_0, nBS_1, nBS_2 )) = user_per_tier(1+getTier( uplink_BS_associated(i), nBS_0, nBS_1, nBS_2 )) + 1;
end

uplink_nU_BS_0 = user_per_tier(1);
uplink_nU_BS_1 = user_per_tier(2);
uplink_nU_BS_2 = user_per_tier(3);

end