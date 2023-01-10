function [Pr_dB, power_dB] = getNewReceivedPower(Pt_dB, association_matrix, U, nU, BSLocation, BSType, nBS_0, nBS_1, nBS_2, channel_dB)

Ux = (U(:, 1));
Uy = (U(:, 2));
BSx = BSLocation(:, 1);
BSy = BSLocation(:, 2);

Pr_dB = zeros(nU, length(BSx));
L_dB = zeros(nU, length(BSx));

% for b = 1 : nBS_0 + nBS_1 + nBS_2
%     users_associated = find(association_matrix(:, b)==1);
%     d(:, b) = sqrt( (Ux - BSx(b)).^2 + (Uy - BSy(b)).^2 );
%     L_dB(:, b) = getPathLoss(BSType(b), d(:, b));
%     Pr_dB(:, b) = Pt_dB(:, b) - L_dB(:, b) + channel_dB(:, b);
%     power_dB(b) = Pr_dB(users_associated, b);
% end

% for i = 1 : nU
%     d(i, :) = sqrt( (Ux(i) - BSx).^2 + (Uy(i) - BSy).^2 );
%     L_dB(i, :) = getPathLoss(BSType, d(i, :));
%     Pr_dB(i, :) = Pt_dB(1, :) - (L_dB(i, :)) + channel_dB(i, :);
%     power_dB(i) = Pr_dB(i, BS_associated(i));    
% end





for i = 1 : length(U)
   BS_associated(i) = find(association_matrix(i, :)==1);
    d(i, :) = sqrt( (Ux(i) - BSx).^2 + (Uy(i) - BSy).^2 );
    L_dB(i, :) = getPathLoss(BSType, d(i, :));
    Pr_dB(i, :) = Pt_dB(1, :) - L_dB(i, :) + channel_dB(i, :);
    power_dB(i) = Pr_dB(i, BS_associated(i));
end



end

