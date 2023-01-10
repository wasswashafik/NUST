clc;
clear all
close all
% -------------------------------------------------------------------------

area = 1000*1000;                  % the total area in squared meters

user_intensity = 1e-4;          % the intensity of users in m^2
UHFM_intensity = 5e-6;              % the BS intensity in m^2
ratio = 5;
mm_intensity = ratio*UHFM_intensity;
%deploy users in the area

SINR_Threshold = -10;
R_min_downlink = 0:0.1e6:20e6; % minimum rate is 0.5 Mbps
R_min_uplink = 0:0.1e6:20e6;
Interference_threshold_downlink = 10^-7;
Interference_threshold_dB_downlink = 10*log10(Interference_threshold_downlink);
Interference_threshold_uplink = 10^-10;
Interference_threshold_dB_uplink = 10*log10(Interference_threshold_uplink);

Niterations = 500;
warning off;
beam_width = 10;    %beam width in degrees

plotable_sum_rate_downlink = zeros([1 length(R_min_downlink)]);
plotable_uhf_macro_average_rate_per_user_downlink = zeros([1 length(R_min_downlink)]);
plotable_mm_small_average_rate_per_user_downlink = zeros([1 length(R_min_downlink)]);
plotable_uhf_small_average_rate_per_user_downlink = zeros([1 length(R_min_downlink)]);
plotable_system_coverage_downlink = zeros([1 length(R_min_downlink)]);
plotable_tier0_coverage_downlink = zeros([1 length(R_min_downlink)]);
plotable_tier1_coverage_downlink = zeros([1 length(R_min_downlink)]);
plotable_tier2_coverage_downlink = zeros([1 length(R_min_downlink)]);

plotable_sum_rate_uplink = zeros([1 length(R_min_uplink)]);
plotable_uhf_macro_average_rate_per_user_uplink = zeros([1 length(R_min_downlink)]);
plotable_mm_small_average_rate_per_user_uplink = zeros([1 length(R_min_downlink)]);
plotable_uhf_small_average_rate_per_user_uplink = zeros([1 length(R_min_downlink)]);
plotable_system_coverage_uplink = zeros([1 length(R_min_downlink)]);
plotable_tier0_coverage_uplink = zeros([1 length(R_min_downlink)]);
plotable_tier1_coverage_uplink = zeros([1 length(R_min_downlink)]);
plotable_tier2_coverage_uplink = zeros([1 length(R_min_downlink)]);

for p = 1 : length(R_min_downlink)
    
    temp_sumrate_down = 0;      temp_sumrate_up = 0;
    temp_syscov_down = 0;       temp_syscov_up = 0;
    temp_0_cov_down = 0;        temp_0_cov_up = 0;
    temp_1_cov_down = 0;        temp_1_cov_up = 0;
    temp_2_cov_down = 0;        temp_2_cov_up = 0;
    
    temp_uhf_macro_average_rate_per_user_downlink = 0;
    temp_mm_small_average_rate_per_user_downlink = 0;
    temp_uhf_small_average_rate_per_user_downlink = 0;
    
    temp_uhf_macro_average_rate_per_user_uplink = 0;
    temp_mm_small_average_rate_per_user_uplink = 0;
    temp_uhf_small_average_rate_per_user_uplink = 0;
    
    
    for j = 1:Niterations
        
        [nBS_0, nBS_1, nBS_2, BSType, BSLocation] = deployBS(area, UHFM_intensity, ratio, 0.1);
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
        temp_sum_rate_per_tier_downlink = zeros([1 3]);
        rate_downlink = zeros([1 nU]);
        
        [Pt_dB_downlink, B_downlink, f_downlink, subcarriers_downlink] = getBSProperties(BSType);
        [nU_BS_0_downlink, nU_BS_1_downlink, nU_BS_2_downlink, Pr_dB_downlink, power_dB_downlink, L_dB_downlink, d, association_matrix_downlink, tier_per_user_downlink, BS_associated_downlink] = downlinkUserAssociation(U, nU, Pt_dB_downlink, nBS_0, nBS_1, nBS_2, BSType, BSLocation);
        [Pt_dB_downlink_new] = downlinkPowerControl(U, nU, nBS_0, nBS_1, nBS_2, BSType, BSLocation, association_matrix_downlink, Pt_dB_downlink, Interference_threshold_dB_downlink, channel_dB);
        [Pr_dB_downlink, power_dB_downlink] = getNewReceivedPower(Pt_dB_downlink_new, association_matrix_downlink, U, nU, BSLocation, BSType, nBS_0, nBS_1, nBS_2, channel_dB);
        bandwidth_array_downlink = getDownlinkBandwidthPerUser(association_matrix_downlink, BSType, BS_associated_downlink);
        [SINR_downlink] = downlinkSINR(bandwidth_array_downlink, nU, power_dB_downlink, Pr_dB_downlink, Pt_dB_downlink_new, BS_associated_downlink, BSType, nBS_0, nBS_1, nBS_2, association_matrix_downlink, angle_array, beam_width );
        rate_downlink = bandwidth_array_downlink.* log2(1+(10.^(SINR_downlink/10))).*coverage_matrix_downlink;
        
        for f = 1 : nU
            if rate_downlink(f) < R_min_downlink(p)
                coverage_matrix_downlink(f) = 0;
            else
                covered_users_per_tier_downlink(1+BSType(BS_associated_downlink(f))) = covered_users_per_tier_downlink(1+BSType(BS_associated_downlink(f)))+1;
                temp_sum_rate_per_tier_downlink(1+BSType(BS_associated_downlink(f))) = temp_sum_rate_per_tier_downlink(1+BSType(BS_associated_downlink(f))) + rate_downlink(f);
            end
        end
        
        coverage_probability = sum(coverage_matrix_downlink/nU);
        system_rate = sum(rate_downlink.*coverage_matrix_downlink);
        
        temp_sumrate_down = temp_sumrate_down + system_rate;
        temp_uhf_macro_average_rate_per_user_downlink = temp_uhf_macro_average_rate_per_user_downlink + temp_sum_rate_per_tier_downlink(1)./covered_users_per_tier_downlink(1);
        temp_mm_small_average_rate_per_user_downlink = temp_mm_small_average_rate_per_user_downlink + temp_sum_rate_per_tier_downlink(2)./covered_users_per_tier_downlink(2);
        temp_uhf_small_average_rate_per_user_downlink = temp_uhf_small_average_rate_per_user_downlink + temp_sum_rate_per_tier_downlink(3)./covered_users_per_tier_downlink(3);
        
        temp_syscov_down = temp_syscov_down + coverage_probability ;
        temp_0_cov_down = temp_0_cov_down + (covered_users_per_tier_downlink(1)/nU) ;
        temp_1_cov_down = temp_1_cov_down + (covered_users_per_tier_downlink(2)/nU) ;
        temp_2_cov_down = temp_2_cov_down + (covered_users_per_tier_downlink(3)/nU) ;
        
        
        % --------------------------------------------------------------------------------------------------------------
        % --------------------------------------------------------------------------------------------------------------
        %                                   UPLINK
        % --------------------------------------------------------------------------------------------------------------
        % --------------------------------------------------------------------------------------------------------------
        
        covered_users_per_tier_uplink = zeros([1 3]);
        temp_sum_rate_per_tier_uplink = zeros([1 3]);
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
        rate_uplink =  bandwidth_array_downlink.* log2(1+(10.^(SINR_uplink/10))).*coverage_matrix_uplink;
        for Us = 1:nU
            if rate_uplink(Us) < R_min_uplink(p)
                coverage_matrix_uplink(Us) = 0;
            else
                covered_users_per_tier_uplink(1+BSType(BS_associated_downlink(Us))) = covered_users_per_tier_uplink(1+BSType(BS_associated_downlink(Us))) + 1;
                temp_sum_rate_per_tier_uplink(1+BSType(BS_associated_downlink(Us))) = temp_sum_rate_per_tier_uplink(1+BSType(BS_associated_downlink(f))) + rate_uplink(Us);
            end
        end
        
        coverage_probability_uplink = sum(coverage_matrix_uplink/nU);
        system_rate_uplink = sum(rate_uplink.*coverage_matrix_uplink);
        
        temp_uhf_macro_average_rate_per_user_uplink = temp_uhf_macro_average_rate_per_user_uplink + (temp_sum_rate_per_tier_uplink(1))./(covered_users_per_tier_uplink(1));
        temp_mm_small_average_rate_per_user_uplink = temp_uhf_macro_average_rate_per_user_uplink + (temp_sum_rate_per_tier_uplink(2))./(covered_users_per_tier_uplink(2));
        temp_uhf_small_average_rate_per_user_uplink = temp_uhf_macro_average_rate_per_user_uplink + (temp_sum_rate_per_tier_uplink(3))./(covered_users_per_tier_uplink(3));
        
        temp_sumrate_up = temp_sumrate_up + system_rate_uplink;
        temp_syscov_up = temp_syscov_up + coverage_probability_uplink ;
        temp_0_cov_up = temp_0_cov_up + (covered_users_per_tier_uplink(1)/nU) ;
        temp_1_cov_up = temp_1_cov_up + (covered_users_per_tier_uplink(2)/nU) ;
        temp_2_cov_up = temp_2_cov_up + (covered_users_per_tier_uplink(3)/nU) ;
        
        
        clc;
        display(strcat('Percent Completion:  ', strcat(num2str(ceil(100*p/length(ratio)))),'%'));
        display(strcat('Percent Completion:  ', strcat(num2str(ceil(100*j/Niterations))),'%'));
        
    end
    
    plotable_uhf_macro_average_rate_per_user_downlink(p) = temp_uhf_macro_average_rate_per_user_downlink/Niterations;
    plotable_mm_small_average_rate_per_user_downlink(p) = temp_mm_small_average_rate_per_user_downlink/Niterations;
    plotable_uhf_small_average_rate_per_user_downlink(p) = temp_uhf_small_average_rate_per_user_downlink/Niterations;
    plotable_sum_rate_downlink(p) = temp_sumrate_down/Niterations;
    plotable_system_coverage_downlink(p) = temp_syscov_down/Niterations;
    plotable_tier0_coverage_downlink(p) = temp_0_cov_down/Niterations;
    plotable_tier1_coverage_downlink(p) = temp_1_cov_down/Niterations;
    plotable_tier2_coverage_downlink(p) = temp_2_cov_down/Niterations;
    
    
    plotable_uhf_macro_average_rate_per_user_uplink(p) = temp_uhf_macro_average_rate_per_user_uplink/Niterations;
    plotable_mm_small_average_rate_per_user_uplink(p) = temp_mm_small_average_rate_per_user_uplink/Niterations;
    plotable_uhf_small_average_rate_per_user_uplink(p) = temp_uhf_small_average_rate_per_user_uplink/Niterations;
    plotable_sum_rate_uplink(p) = temp_sumrate_up/Niterations;
    plotable_system_coverage_uplink(p) = temp_syscov_up/Niterations;
    plotable_tier0_coverage_uplink(p) = temp_0_cov_up/Niterations;
    plotable_tier1_coverage_uplink(p) = temp_1_cov_up/Niterations;
    plotable_tier2_coverage_uplink(p) = temp_2_cov_up/Niterations;
    
    coupled_users_fraction(p) = sum(sum(association_matrix_downlink.*association_matrix_uplink))/nU;
    decoupled_users_fraction(p) = 1 - coupled_users_fraction(p);
end

figure;
plot(R_min_downlink, plotable_tier0_coverage_downlink, 'r');
title('Coverage Probability - Downlink ');
xlabel('R_m_i_n (bps) -->');
ylabel('Coverage Probability');
hold on;
plot(R_min_downlink, plotable_tier1_coverage_downlink, 'g');
plot(R_min_downlink, plotable_tier2_coverage_downlink, 'k');
legend('UHF Macro', 'mmWave Small', 'UHF Small');

figure; 
plot(R_min_downlink, plotable_system_coverage_downlink, 'm');
title('Coverage Probability - Downlink ');
xlabel('R_m_i_n (bps) -->');
ylabel('Coverage Probability');

% figure;
% plot(R_min_downlink, plotable_sum_rate_downlink, 'r');
% title('Rate - Downlink');
% xlabel('l_s / l_u_h_f');
% ylabel('Rate(bps)');
% figure; 
% title('Rate per user - Downlink');
% xlabel('l_s / l_u_h_f');
% ylabel('Rate(bps)');
% plot(R_min_downlink, plotable_uhf_macro_average_rate_per_user_downlink, 'g');
% hold on
% plot(R_min_downlink, plotable_mm_small_average_rate_per_user_downlink, 'k');
% plot(R_min_downlink, plotable_uhf_small_average_rate_per_user_downlink, 'm');
% legend('UHF Macro average rate per user', 'mmWave Small average rate per user', 'UHF Small average rate per user');




figure;
plot(R_min_uplink, plotable_tier0_coverage_uplink, 'r');
title('Coverage Probability - Uplink ');
xlabel('R_m_i_n (bps) -->');
ylabel('Coverage Probability');
hold on
plot(R_min_uplink, plotable_tier1_coverage_uplink, 'g');
hold on;
plot(R_min_uplink, plotable_tier2_coverage_uplink, 'k');
legend('UHF Macro', 'mmWave Small', 'UHF Small');

figure;
plot(R_min_uplink, plotable_system_coverage_uplink, 'm');
title('Coverage Probability - Uplink ');
xlabel('R_m_i_n (bps) -->');
ylabel('Coverage Probability');



% figure;
% plot(R_min_uplink, plotable_sum_rate_uplink, 'r');
% title('Rate - Uplink');
% xlabel('l_s / l_u_h_f');
% ylabel('Rate(bps)');
% figure;
% title('Rate per user - Uplink');
% xlabel('l_s / l_u_h_f');
% ylabel('Rate(bps)');
% plot(R_min_uplink, plotable_uhf_macro_average_rate_per_user_uplink, 'g');
% hold on
% plot(R_min_uplink, plotable_mm_small_average_rate_per_user_uplink, 'k');
% plot(R_min_uplink, plotable_uhf_small_average_rate_per_user_uplink, 'm');
% legend('UHF Macro average rate per user', 'mmWave Small average rate per user', 'UHF Small average rate per user');


% figure;
% plot(R_min_uplink, coupled_users_fraction, ratio, decoupled_users_fraction);
% title('Decoupling');
% xlabel('l_s / l_u_h_f');
% ylabel('Fraction of users');
% legend('Coupled', 'Decoupled');
% display('Done');



