function [nU_BS_0,nU_BS_1,nU_BS_2,Pr_dB, power_dB, L_dB, d, association_matrix, downlink_tier_per_user, BS_associated] = downlinkUserAssociation(U, nU, Pt_dB, nBS_0, nBS_1, nBS_2, BSType, BSLocation)
%These lines are redundant and can be removed
Ux = (U(:, 1));
Uy = (U(:, 2));
BSx = BSLocation(:, 1);
BSy = BSLocation(:, 2);
%--------------------------------------------
nU_BS_0 = 0;
nU_BS_1 = 0;
nU_BS_2 = 0;
%Defining and intializing needed variables
Pr_dB = zeros(nU, length(BSx));
L_dB = zeros(nU, length(BSx));
power_dB = zeros(nU,1);
BS_associated = zeros(nU,1);
association_matrix = zeros(nU, length(BSx));
d = zeros(nU, length(BSx));
tier = zeros(nU,1);

%Adding antenna gain for mmWave cells
antenna_gain = [getDownlinkAntennaGain(0, 1).*ones([1 nBS_0])  getDownlinkAntennaGain(1, 1).*ones([1 nBS_1])   getDownlinkAntennaGain(2, 1).*ones([1 nBS_2])];  %add antenna gain

for i = 1 : nU
    %Calculate distance of each user and BS, get distance array
    d(i, :) = sqrt( (Ux(i) - BSx).^2 + (Uy(i) - BSy).^2 );
    
    %Use distance to find pathloss
    L_dB(i, :) = getPathLoss(BSType, d(i, :));
    
    %Get power recieved by subtracting pathloss
    Pr_dB(i, :) = Pt_dB(1, :) - L_dB(i, :);
    
    %Add antenna gain to mmWave cells only, this adds 18dB
    Pr_dB(i, :) = Pr_dB(i, :) + antenna_gain;
    
    %Add bias to mmWave cells, for that get bias array
    bias_array = [getBiasFactor(0).*ones([1 nBS_0])  getBiasFactor(1).*ones([1 nBS_1])   getBiasFactor(2).*ones([1 nBS_2])];

    %Find BS with maximum power and get its power and index to associate it
    %with the user
    [power_dB(i), BS_associated(i)] = max(Pr_dB(i, :)+bias_array);
    
    %Remove the bias
    power_dB(i) = power_dB(i) - getBiasFactor(getTier(BS_associated(i),nBS_0, nBS_1, nBS_2));
    
    %create association matrix
    association_matrix(i, BS_associated(i)) = 1;
    
    %get tier of the BS connected to each user
    tier(i) = getTier(BS_associated(i), nBS_0, nBS_1, nBS_2);
    
    
    if tier(i) == 0
        nU_BS_0 = nU_BS_0 + 1;
    else if tier(i) == 1
            nU_BS_1 = nU_BS_1 + 1;
        else
            nU_BS_2 = nU_BS_2 + 1;
        end
    end
    
    
    downlink_tier_per_user = tier;
end

end