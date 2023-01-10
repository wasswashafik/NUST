function [ P ] = downlinkPowerControl(U, nU, nBS_0, nBS_1, nBS_2, BSType, BSLocation, association_matrix_downlink, Pt_dB_downlink, Interference_threshold_dB, channel_dB)

P = Pt_dB_downlink;
It = Interference_threshold_dB;
Ux = (U(:, 1));
Uy = (U(:, 2));
BSx = BSLocation(:, 1);
BSy = BSLocation(:, 2);
d = zeros(nU, length(BSLocation));
BS_associated = zeros(1, nU);

for u = 1 : length(U)
    p = association_matrix_downlink;
    BS_associated(u) = find(p(u, :) == 1);
    d(u, :) = sqrt((Ux(u) - BSx).^2 + (Uy(u) - BSy).^2);
    if BSType(BS_associated(u)) == 0
        p = [association_matrix_downlink(:, 1:nBS_0), association_matrix_downlink(:, (nBS_0+nBS_1+1):(nBS_0+nBS_1+nBS_2))];
    else if BSType(BS_associated(u)) == 1
            continue
            %p = association_matrix_downlink(:, nBS_0+1:nBS_0+nBS_1);
        else if BSType(BS_associated(u)) == 2
                p = [association_matrix_downlink(:, 1:nBS_0), association_matrix_downlink(:, (nBS_0+nBS_1+1):length(BSLocation))];
            end
        end
    end
    
    l = zeros(nU, size(p, 2)-1);
    if isempty(find(p(u, :) == 0))
        %error('No interfering BS');
        continue;
    end
    l(u, :) = transpose(sqrt( (Ux(u)-BSx(find(p(u, :) == 0))).^2 + (Uy(u)-BSy(find(p(u, :) == 0))).^2 ));
%    L_dB(u, :) = getPathLoss(BSType, l(u, :));
    min_distance(u) = min(l(u, :));
    
    nearest_interfering_BS(u) = find(d(u, :) == min_distance(u));
    PL(u) = getPathLoss(BSType(nearest_interfering_BS(u)), min_distance(u));
    P(nearest_interfering_BS(u)) = min(P(nearest_interfering_BS(u)), It + PL(u) - channel_dB(u, nearest_interfering_BS(u)));
end
end

