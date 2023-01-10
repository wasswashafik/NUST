function [ L_dB ] = getPathLoss( tier, distance, blocking_prob)

%using pathloss formula, distance, and tier, calculate pathloss

alpha_tier0_1 = 2; alpha_tier0_2 = 3;
alpha_tier1_LoS_1 = 2; alpha_tier1_LoS_2 = 2.2;
alpha_tier1_NLoS_1 = 3.3; alpha_tier1_NLoS_2 = 3.5;
alpha_tier2_1 = 2; alpha_tier2_2 = 3;
critical_radius_tier0 = 500;%375
critical_radius_tier1 = 100;
critical_radius_tier2 = 200;

L_dB = zeros(length(tier), 1);
c = 3e8;

%Get pathloss exponent
%PLE = assignPLE(tier, distance,blocking_prob);
%refernece distance
%redundant variable
do = 1;

%path loss =           fixed                 +         distance dependent


%PLE = assignPLE(tier, distance,blocking_prob);
%L_dB = 20.*log10(4*pi*getFrequency(tier)./c) + (10.*PLE).*log10((distance)./do);

for k = 1 : length(tier)
    if tier(k) == 0
        if distance(k) < critical_radius_tier0
            L_dB(k) = 20.*log10(4*pi*getFrequency(tier(k))./c) + (10.*alpha_tier0_1).*log10((distance(1, k))./do);
        else
            L_dB(k) = 20.*log10(4*pi*getFrequency(tier(k))./c) + (10.*alpha_tier0_1).*log10((critical_radius_tier0)./do) + (10.*alpha_tier0_2).*log10((distance(1, k))./critical_radius_tier0);
        end
    else if tier(k) == 1
            if distance < critical_radius_tier1
                if isLOS(distance(k),blocking_prob)
                    L_dB(k) = 20.*log10(4*pi*getFrequency(tier(k))./c) + (10.*alpha_tier1_LoS_1).*log10((distance(1, k))./do);
                else
                    L_dB(k) = 20.*log10(4*pi*getFrequency(tier(k))./c) + (10.*alpha_tier1_NLoS_1).*log10((distance(1, k))./do);
                end
            else
                if isLOS(distance(k),blocking_prob)
                    L_dB(k) = 20.*log10(4*pi*getFrequency(tier(k))./c) + (10.*alpha_tier1_LoS_1).*log10((critical_radius_tier1)./do) + (10.*alpha_tier1_LoS_2).*log10(distance(1, k)./critical_radius_tier1);
                else
                    L_dB(k) = 20.*log10(4*pi*getFrequency(tier(k))./c) + (10.*alpha_tier1_NLoS_2).*log10((critical_radius_tier1)./do) + (10.*alpha_tier1_NLoS_2).*log10((distance(1, k))./critical_radius_tier1);
                end
            end
        else if tier(k) == 2
                if distance(k) < critical_radius_tier2
                    L_dB(k) = 20.*log10(4*pi*getFrequency(tier(k))./c) + (10.*alpha_tier2_1).*log10((distance(1, k))./do);
                else
                    L_dB(k) = 20.*log10(4*pi*getFrequency(tier(k))./c) + (10.*alpha_tier2_1).*log10((critical_radius_tier2)./do) + (10.*alpha_tier0_2).*log10((distance(1, k))./critical_radius_tier2);
                end
            end
        end
    end
end
end

