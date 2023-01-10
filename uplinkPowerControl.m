function [P] = uplinkPowerControl(U, nU, nBS_0, nBS_1, nBS_2, BSType, BSLocation, association_matrix_uplink, Pt_dB_uplink, Interference_threshold_dB_uplink, channel_dB)

P = Pt_dB_uplink;
It = Interference_threshold_dB_uplink;
Ux = (U(:, 1));
Uy = (U(:, 2));
BSx = BSLocation(:, 1);
BSy = BSLocation(:, 2);
d = zeros(nU, nBS_0+nBS_1+nBS_2);
PL = zeros(1, nBS_0+nBS_1+nBS_2);
for b = 1 : nBS_0+nBS_1+nBS_2
    if ((b >=nBS_0+1) && (b<=nBS_0+nBS_1))
        continue
    end
    p = association_matrix_uplink;
    %users_associated(b, :) = find(p(:, b)==1);
    d(:, b) = sqrt((Ux - BSx(b)).^2 + (Uy - BSy(b)).^2);
    if isempty(find(p(:, b) == 0))
        %error('No interfering BS');
        continue;
    end
    %l = zeros(length(find(p(:, b)==1)), length(BSLocation));
    %l(:, b) = transpose(sqrt( (Ux(find(p(:, b)==0))-BSx(b)).^2 + (Uy(find(p(:, b)==0))-BSy(b)).^2 ));
    l = transpose(sqrt( (Ux(find(p(:, b)==0))-BSx(b)).^2 + (Uy(find(p(:, b)==0))-BSy(b)).^2 ));
    min_distance = min(l);
    
% %     if BSType(b) == 0
% %         p = [association_matrix_uplink(:, 1:nBS_0), association_matrix_uplink(:, (nBS_0+nBS_1+1):(nBS_0+nBS_1+nBS_2))];
% %     else if BSType(b) == 1
% %             p = association_matrix_uplink(:, nBS_0+1:nBS_0+nBS_1);
% %         else if BSType(b) == 2
% %                 p = [association_matrix_uplink(:, 1:nBS_0), association_matrix_uplink(:, (nBS_0+nBS_1+1):(nBS_0+nBS_1+nBS_2))];
% %             end
% %         end
% %     end
    nearest_interfering_user(b) = find(d(:, b) == min_distance);
    PL(b) = getPathLoss(BSType(b), min_distance);
    P(nearest_interfering_user(b)) = min(P(nearest_interfering_user(b)), It + PL(b) - channel_dB(nearest_interfering_user(b), b));
        
end


end