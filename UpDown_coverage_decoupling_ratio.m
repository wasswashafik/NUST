clc;
clear all
close all
% -------------------------------------------------------------------------

general_filename = 'X:\FYP\TEST CODE\NUST_small_to_large\Data';


area = 1000*1000;                  % the total area in squared meters
%area = 500*100;
user_intensity = 1e-4;          % the intensity of users in m^2
UHFM_intensity = 5e-6;              % the BS intensity in m^2
ratio =  0:1:10;
mm_intensity = ratio*UHFM_intensity;
%deploy users in the area

%SINR_Threshold = -10;
R_min_downlink = 0.5e6; % minimum rate is 1 Mbps
R_min_uplink = 0.5e5; % minimum threshold rate for
Interference_threshold_downlink = 10^-10;
Interference_threshold_dB_downlink = 10*log10(Interference_threshold_downlink);
Interference_threshold_uplink = 10^-10;
Interference_threshold_dB_uplink = 10*log10(Interference_threshold_uplink);

Niterations = 10;
warning off;
beam_width = 10;    % beam width in degrees

plotable_sum_rate_downlink = zeros([1 length(ratio)]);
plotable_uhf_macro_average_rate_per_user_downlink = zeros([1 length(ratio)]);
plotable_mm_small_average_rate_per_user_downlink =  zeros([1 length(ratio)]);
plotable_uhf_small_average_rate_per_user_downlink = zeros([1 length(ratio)]);
plotable_system_coverage_downlink = zeros([1 length(ratio)]);
plotable_tier0_coverage_downlink =  zeros([1 length(ratio)]);
plotable_tier1_coverage_downlink =  zeros([1 length(ratio)]);
plotable_tier2_coverage_downlink =  zeros([1 length(ratio)]);
plotable_tier0_sum_rate_downlink = zeros([1 length(ratio)]);
plotable_tier1_sum_rate_downlink = zeros([1 length(ratio)]);
plotable_tier2_sum_rate_downlink = zeros([1 length(ratio)]);

plotable_sum_rate_uplink = zeros([1 length(ratio)]);
plotable_uhf_macro_average_rate_per_user_uplink = zeros([1 length(ratio)]);
plotable_mm_small_average_rate_per_user_uplink =  zeros([1 length(ratio)]);
plotable_uhf_small_average_rate_per_user_uplink = zeros([1 length(ratio)]);
plotable_system_coverage_uplink = zeros([1 length(ratio)]);
plotable_tier0_sum_rate_uplink = zeros([1 length(ratio)]);
plotable_tier1_sum_rate_uplink = zeros([1 length(ratio)]);
plotable_tier2_sum_rate_uplink = zeros([1 length(ratio)]);
plotable_tier0_coverage_uplink =  zeros([1 length(ratio)]);
plotable_tier1_coverage_uplink =  zeros([1 length(ratio)]);
plotable_tier2_coverage_uplink =  zeros([1 length(ratio)]);

plotable_system_coverage = zeros([1 length(ratio)]);
plotable_coupled_users_fraction = zeros([1 length(ratio)]);
plotable_decoupled_users_fraction = zeros([1 length(ratio)]);

for p = 1 : length(ratio)
    
    temp_sumrate_down = 0;      temp_sumrate_up = 0;
    temp_syscov_down  = 0;       temp_syscov_up = 0;
    temp_0_cov_down   = 0;        temp_0_cov_up = 0;
    temp_1_cov_down   = 0;        temp_1_cov_up = 0;
    temp_2_cov_down   = 0;        temp_2_cov_up = 0;
    temp_system_coverage = 0;
    
    temp_uhf_macro_average_rate_per_user_downlink = 0;
    temp_mm_small_average_rate_per_user_downlink  = 0;
    temp_uhf_small_average_rate_per_user_downlink = 0;
    temp_sum_rate_per_tier_downlink = zeros([1 3]);
    
    temp_sum_rate_per_tier_uplink = zeros([1 3]);
    temp_uhf_macro_average_rate_per_user_uplink = 0;
    temp_mm_small_average_rate_per_user_uplink  = 0;
    temp_uhf_small_average_rate_per_user_uplink = 0;
    
    temp_coupled_users_fraction   = 0;
    temp_decoupled_users_fraction = 0;
    
    for j = 1:Niterations
        
        [nBS_0, nBS_1, nBS_2, BSType, BSLocation] = deployBS(area, UHFM_intensity, ratio(p), 0.2);
        [U, nU] = deployUsers(user_intensity, area);
        angle_array = getBsToUserAngle(U, nU, nBS_0, nBS_1, BSLocation);
        [channel, channel_dB] = getChannel(nU, BSType);
        
        % --------------------------------------------------------------------------------------------------------------
        % --------------------------------------------------------------------------------------------------------------
        %                                   DOWNLINK
        % --------------------------------------------------------------------------------------------------------------
        % --------------------------------------------------------------------------------------------------------------
        
        coverage_probability_downlink = 0;
        coverage_matrix_downlink = ones([1 nU]);
        covered_users_per_tier_downlink = zeros([1 3]);
        rate_downlink = zeros([1 nU]);
        
        [Pt_dB_downlink, B_downlink, f_downlink, subcarriers_downlink] = getBSProperties(BSType);
        [nU_BS_0_downlink, nU_BS_1_downlink, nU_BS_2_downlink, Pr_dB_downlink, power_dB_downlink, L_dB_downlink, d, association_matrix_downlink, tier_per_user_downlink, BS_associated_downlink] = downlinkUserAssociation(U, nU, Pt_dB_downlink, nBS_0, nBS_1, nBS_2, BSType, BSLocation);
        [Pt_dB_downlink_new] = downlinkPowerControl(U, nU, nBS_0, nBS_1, nBS_2, BSType, BSLocation, association_matrix_downlink, Pt_dB_downlink, Interference_threshold_dB_downlink, channel_dB);
        [Pr_dB_downlink, power_dB_downlink] = getNewReceivedPower(Pt_dB_downlink_new, association_matrix_downlink, U, nU, BSLocation, BSType, nBS_0, nBS_1, nBS_2, channel_dB);
        bandwidth_array_downlink = getDownlinkBandwidthPerUser(association_matrix_downlink, BSType, BS_associated_downlink);
        [SINR_downlink] = downlinkSINR(bandwidth_array_downlink, nU, power_dB_downlink, Pr_dB_downlink, Pt_dB_downlink_new, BS_associated_downlink, BSType, nBS_0, nBS_1, nBS_2, association_matrix_downlink, angle_array, beam_width );
        rate_downlink = bandwidth_array_downlink.* log2(1+(10.^(SINR_downlink/10))).*coverage_matrix_downlink;
        
        
        
        for f = 1 : nU
            if rate_downlink(f) < R_min_downlink
                coverage_matrix_downlink(f) = 0;
            else
                covered_users_per_tier_downlink(1+BSType(BS_associated_downlink(f))) = covered_users_per_tier_downlink(1+BSType(BS_associated_downlink(f)))+1;
                temp_sum_rate_per_tier_downlink(1+BSType(BS_associated_downlink(f))) = temp_sum_rate_per_tier_downlink(1+BSType(BS_associated_downlink(f))) + rate_downlink(f);
            end
        end
        
        
        coverage_probability_downlink = sum(coverage_matrix_downlink/nU);
        system_rate = sum(rate_downlink.*coverage_matrix_downlink);
        temp_sumrate_down = temp_sumrate_down + system_rate;
        
        
        if (covered_users_per_tier_downlink(1) ~= 0)
            temp_uhf_macro_average_rate_per_user_downlink = temp_uhf_macro_average_rate_per_user_downlink + temp_sum_rate_per_tier_downlink(1)./covered_users_per_tier_downlink(1);
        end
        if (covered_users_per_tier_downlink(2) ~= 0)
            temp_mm_small_average_rate_per_user_downlink = temp_mm_small_average_rate_per_user_downlink + temp_sum_rate_per_tier_downlink(2)./covered_users_per_tier_downlink(2);
        end
        if (covered_users_per_tier_downlink(3) ~= 0)
            temp_uhf_small_average_rate_per_user_downlink = temp_uhf_small_average_rate_per_user_downlink + temp_sum_rate_per_tier_downlink(3)./covered_users_per_tier_downlink(3);
        end
        
        
        temp_syscov_down = temp_syscov_down + coverage_probability_downlink ;
        temp_0_cov_down = temp_0_cov_down + (covered_users_per_tier_downlink(1)/nU) ;
        temp_1_cov_down = temp_1_cov_down + (covered_users_per_tier_downlink(2)/nU) ;
        temp_2_cov_down = temp_2_cov_down + (covered_users_per_tier_downlink(3)/nU) ;
        
        
        % --------------------------------------------------------------------------------------------------------------
        % --------------------------------------------------------------------------------------------------------------
        %                                   UPLINK
        % --------------------------------------------------------------------------------------------------------------
        % --------------------------------------------------------------------------------------------------------------
        
        covered_users_per_tier_uplink = zeros([1 3]);
        SINR_uplink = zeros([1 nU]);
        rate_uplink = zeros([1 nU]);
        coverage_matrix_uplink = ones([1 nU]);
        covered_users_uplink = 0;
        coverage_probability_uplink = 0;
        P_max_uplink = 23-30;
        Pt_dB_uplink = P_max_uplink*ones(1, nU);
        [uplink_BS_associated, nU_BS_0_uplink, nU_BS_1_uplink, nU_BS_2_uplink, association_matrix_uplink,  Pr_dB_uplink, max_power_db ] = uplinkUserAssociation( Pt_dB_uplink, U, nU, nBS_0, nBS_1, nBS_2, BSType, BSLocation);
        [Pt_dB_uplink_new] = uplinkPowerControl(U, nU, nBS_0, nBS_1, nBS_2, BSType, BSLocation, association_matrix_uplink, Pt_dB_uplink, Interference_threshold_dB_uplink, channel_dB);
        [Pr_dB_uplink_new, power_dB_uplink] = getNewReceivedPower(transpose(Pt_dB_uplink_new), association_matrix_uplink, U, nU, BSLocation, BSType, nBS_0, nBS_1, nBS_2, channel_dB);
        bandwidth_array_uplink = getDownlinkBandwidthPerUser(association_matrix_uplink, BSType, uplink_BS_associated);
        [uplink_interference_watts, uplink_interference_dB, IplusN_db] = getUplinkInterference( uplink_BS_associated,bandwidth_array_uplink,association_matrix_uplink,BSType, Pr_dB_uplink_new, power_dB_uplink, nU);
        SINR_uplink = (power_dB_uplink) - IplusN_db;
        rate_uplink =  bandwidth_array_uplink.* log2(1+(10.^(SINR_uplink/10))).*coverage_matrix_uplink;
        
        for Us = 1:nU
            if rate_uplink(Us) < R_min_uplink
                coverage_matrix_uplink(Us) = 0;
            else
                covered_users_per_tier_uplink(1+BSType(BS_associated_downlink(Us))) = covered_users_per_tier_uplink(1+BSType(BS_associated_downlink(Us))) + 1;
                temp_sum_rate_per_tier_uplink(1+BSType(BS_associated_downlink(Us))) = temp_sum_rate_per_tier_uplink(1+BSType(BS_associated_downlink(Us))) + rate_uplink(Us);
            end
        end
        coverage_probability_uplink = sum(coverage_matrix_uplink/nU);
        system_rate_uplink = sum(rate_uplink.*coverage_matrix_uplink);
        
        
        
        if (covered_users_per_tier_uplink(1) ~= 0)
            temp_uhf_macro_average_rate_per_user_uplink = temp_uhf_macro_average_rate_per_user_uplink + (temp_sum_rate_per_tier_uplink(1))./(covered_users_per_tier_uplink(1));
        end
        if (covered_users_per_tier_uplink(2) ~= 0)
            temp_mm_small_average_rate_per_user_uplink = temp_mm_small_average_rate_per_user_uplink + (temp_sum_rate_per_tier_uplink(2))./(covered_users_per_tier_uplink(2));
        end
        if (covered_users_per_tier_uplink(3) ~= 0)
            temp_uhf_small_average_rate_per_user_uplink = temp_uhf_small_average_rate_per_user_uplink + (temp_sum_rate_per_tier_uplink(3))./(covered_users_per_tier_uplink(3));
        end
        temp_sumrate_up = temp_sumrate_up + system_rate_uplink;
        temp_syscov_up = temp_syscov_up + coverage_probability_uplink ;
        temp_0_cov_up = temp_0_cov_up + (covered_users_per_tier_uplink(1)/nU) ;
        temp_1_cov_up = temp_1_cov_up + (covered_users_per_tier_uplink(2)/nU) ;
        temp_2_cov_up = temp_2_cov_up + (covered_users_per_tier_uplink(3)/nU) ;
        temp_system_coverage = temp_system_coverage + (sum(coverage_matrix_downlink.*coverage_matrix_uplink)/nU);
        temp_coupled_users_fraction = temp_coupled_users_fraction + sum(sum(association_matrix_downlink.*association_matrix_uplink))/nU;
        temp_decoupled_users_fraction = temp_decoupled_users_fraction + (1 - sum(sum(association_matrix_downlink.*association_matrix_uplink))/nU);
        
        
        
        clc;
        display(strcat('Percent Completion:  ', strcat(num2str(ceil(100*p/length(ratio)))),'%'));
        display(strcat('Percent Completion:  ', strcat(num2str(ceil(100*j/Niterations))),'%'));
        
    end
    
    
    plotable_uhf_macro_average_rate_per_user_downlink(p) = temp_uhf_macro_average_rate_per_user_downlink/Niterations;
    plotable_mm_small_average_rate_per_user_downlink(p) = temp_mm_small_average_rate_per_user_downlink/Niterations;
    plotable_uhf_small_average_rate_per_user_downlink(p) = temp_uhf_small_average_rate_per_user_downlink/Niterations;
    
    plotable_system_coverage(p) = temp_system_coverage/Niterations;
    plotable_sum_rate_downlink(p) = temp_sumrate_down/Niterations;
    plotable_tier0_sum_rate_downlink(p) = temp_sum_rate_per_tier_downlink(1)/Niterations;
    plotable_tier1_sum_rate_downlink(p) = temp_sum_rate_per_tier_downlink(2)/Niterations;
    plotable_tier2_sum_rate_downlink(p) = temp_sum_rate_per_tier_downlink(3)/Niterations;
    plotable_system_coverage_downlink(p) = temp_syscov_down/Niterations;
    plotable_tier0_coverage_downlink(p) = temp_0_cov_down/Niterations;
    plotable_tier1_coverage_downlink(p) = temp_1_cov_down/Niterations;
    plotable_tier2_coverage_downlink(p) = temp_2_cov_down/Niterations;
    
    
    plotable_uhf_macro_average_rate_per_user_uplink(p) = temp_uhf_macro_average_rate_per_user_uplink/Niterations;
    plotable_mm_small_average_rate_per_user_uplink(p) = temp_mm_small_average_rate_per_user_uplink/Niterations;
    plotable_uhf_small_average_rate_per_user_uplink(p) = temp_uhf_small_average_rate_per_user_uplink/Niterations;
    plotable_sum_rate_uplink(p) = temp_sumrate_up/Niterations;
    plotable_system_coverage_uplink(p) = temp_syscov_up/Niterations;
    plotable_tier0_sum_rate_uplink(p) = temp_sum_rate_per_tier_uplink(1)/Niterations;
    plotable_tier1_sum_rate_uplink(p) = temp_sum_rate_per_tier_uplink(2)/Niterations;
    plotable_tier2_sum_rate_uplink(p) = temp_sum_rate_per_tier_uplink(3)/Niterations;

    plotable_tier0_coverage_uplink(p) = temp_0_cov_up/Niterations;
    plotable_tier1_coverage_uplink(p) = temp_1_cov_up/Niterations;
    plotable_tier2_coverage_uplink(p) = temp_2_cov_up/Niterations;
    
    plotable_coupled_users_fraction(p) = temp_coupled_users_fraction/Niterations;
    plotable_decoupled_users_fraction(p) = temp_decoupled_users_fraction/Niterations;
    
end

figure;
plot(ratio, plotable_tier0_coverage_downlink, 'r');
title('Coverage Probability - Downlink '); xlabel('l_s / l_u_h_f'); ylabel('Coverage Probability');
hold on;
plot(ratio, plotable_tier1_coverage_downlink, 'g');
plot(ratio, plotable_tier2_coverage_downlink, 'k');
plot(ratio, plotable_system_coverage_downlink, 'm');
legend('UHF Macro', 'mmWave Small', 'UHF Small', 'Overall');

figure;
plot(ratio, plotable_sum_rate_downlink, 'r');
title('Rate - Downlink');
xlabel('l_s / l_u_h_f');
ylabel('Rate(bps)');

figure;
title('Sum Rate  - Downlink');
xlabel('l_s / l_u_h_f');
ylabel('Rate(bps)');
title('Rate - Downlink');
plot(ratio, plotable_tier0_sum_rate_downlink, 'g');hold on
plot(ratio, plotable_tier1_sum_rate_downlink, 'k');
plot(ratio, plotable_tier2_sum_rate_downlink, 'm');
legend('Tier 0 Sumrate', 'Tier 1 sumrate', 'Tier 2 sumrate');



figure;
plot(ratio, plotable_tier0_coverage_uplink, 'r');
title('Coverage Probability - Uplink ');
xlabel('l_s / l_u_h_f');
ylabel('Coverage Probability');
hold on;
plot(ratio, plotable_tier1_coverage_uplink, 'g');
hold on;
plot(ratio, plotable_tier2_coverage_uplink, 'k');
hold on;
plot(ratio, plotable_system_coverage_uplink, 'm');
legend('UHF Macro', 'mmWave Small', 'UHF Small', 'Overall');

figure;
plot(ratio, plotable_system_coverage, 'm');
title('System Coverage');
xlabel('l_s / l_u_h_f');
figure;
plot(ratio, plotable_sum_rate_uplink, 'r');
title('Rate - Uplink');
xlabel('l_s / l_u_h_f');
ylabel('Rate(bps)');
figure;
plot(ratio, plotable_tier0_sum_rate_uplink, 'g');
hold on
plot(ratio, plotable_tier1_sum_rate_uplink, 'k');
plot(ratio, plotable_tier2_sum_rate_uplink, 'm');
title('Rate per user - Uplink');
xlabel('l_s / l_u_h_f');
ylabel('Rate(bps)');
legend('Tier 0 sum rate', 'Tier 1 sum rate', 'Tier 2 sum rate');


figure;
plot(ratio, plotable_coupled_users_fraction, ratio, plotable_decoupled_users_fraction);
title('Decoupling');
xlabel('l_s / l_u_h_f');
ylabel('Fraction of users');
legend('Coupled', 'Decoupled');
display('Done');

save(strcat(general_filename, '\ratio.mat'), 'ratio');
% saving downlink variables
save(strcat(general_filename, '\plotable_uhf_macro_average_rate_per_user_downlink.mat'), 'plotable_uhf_macro_average_rate_per_user_downlink');
save(strcat(general_filename, '\plotable_mm_small_average_rate_per_user_downlink.mat'), 'plotable_mm_small_average_rate_per_user_downlink');
save(strcat(general_filename, '\plotable_uhf_small_average_rate_per_user_downlink.mat'), 'plotable_uhf_small_average_rate_per_user_downlink');
save(strcat(general_filename, '\plotable_sum_rate_downlink.mat'), 'plotable_sum_rate_downlink');
save(strcat(general_filename, '\plotable_tier0_sum_rate_downlink.mat'), 'plotable_tier0_sum_rate_downlink');
save(strcat(general_filename, '\plotable_tier1_sum_rate_downlink.mat'), 'plotable_tier1_sum_rate_downlink');
save(strcat(general_filename, '\plotable_tier2_sum_rate_downlink.mat'), 'plotable_tier2_sum_rate_downlink');
save(strcat(general_filename, '\plotable_system_coverage_downlink.mat'), 'plotable_system_coverage_downlink');
save(strcat(general_filename, '\plotable_tier0_coverage_downlink.mat'), 'plotable_tier0_coverage_downlink');
save(strcat(general_filename, '\plotable_tier1_coverage_downlink.mat'), 'plotable_tier1_coverage_downlink');
save(strcat(general_filename, '\plotable_tier2_coverage_downlink.mat'), 'plotable_tier2_coverage_downlink');

save(strcat(general_filename, '\plotable_system_coverage.mat'), 'plotable_system_coverage');
%saving uplink variables

save(strcat(general_filename, '\plotable_uhf_macro_average_rate_per_user_uplink.mat'), 'plotable_uhf_macro_average_rate_per_user_uplink');
save(strcat(general_filename, '\plotable_mm_small_average_rate_per_user_uplink.mat'), 'plotable_mm_small_average_rate_per_user_uplink');
save(strcat(general_filename, '\plotable_uhf_small_average_rate_per_user_uplink.mat'), 'plotable_uhf_small_average_rate_per_user_uplink');
save(strcat(general_filename, '\plotable_sum_rate_uplink.mat'), 'plotable_sum_rate_uplink');
save(strcat(general_filename, '\plotable_tier0_sum_rate_uplink.mat'), 'plotable_tier0_sum_rate_uplink');
save(strcat(general_filename, '\plotable_tier1_sum_rate_uplink.mat'), 'plotable_tier1_sum_rate_uplink');
save(strcat(general_filename, '\plotable_tier2_sum_rate_uplink.mat'), 'plotable_tier2_sum_rate_uplink');
save(strcat(general_filename, '\plotable_system_coverage_uplink.mat'), 'plotable_system_coverage_uplink');
save(strcat(general_filename, '\plotable_tier0_coverage_uplink.mat'), 'plotable_tier0_coverage_uplink');
save(strcat(general_filename, '\plotable_tier1_coverage_uplink.mat'), 'plotable_tier1_coverage_uplink');
save(strcat(general_filename, '\plotable_tier2_coverage_uplink.mat'), 'plotable_tier2_coverage_uplink');

%decoupling variables
save(strcat(general_filename, '\plotable_coupled_users_fraction.mat'), 'plotable_coupled_users_fraction');
save(strcat(general_filename, '\plotable_decoupled_users_fraction.mat'), 'plotable_decoupled_users_fraction');