clc;
clear all
close all
% -------------------------------------------------------------------------
%general_filename = 'X:\FYP\TEST CODE\Chicago_small_to_large\Data';

% x = logspace(-1,2);
% y = exp(x);
% loglog(x,y,'-s')
% grid on

PC = ispc;
UNIX = isunix;
if PC == 1
    seperator = '\';
    
elseif UNIX == 1
    
    seperator = '/';
end

general_filename = strcat(pwd, strcat(seperator, 'Data_NUST'));

mkdir(general_filename)

area = 1000*1000;                  % the total area in squared meters

user_intensity = 0.5e-4;          % the intensity of users in m^2
mmtoUHFFactor = 0.3; % 0 - all mm, 1 - all UHF
UHFM_intensity =  9.5493e-7;              % the BS intensity in m^2
ratio =  10;
mm_intensity = ratio*UHFM_intensity;
%deploy users in the area

SINR_Threshold = 0;
SINR_Threshold_Uplink = -20;
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



plotable_sum_rate_uplink_ND = zeros([1 length(ratio)]);
plotable_uhf_macro_average_rate_per_user_uplink_ND = zeros([1 length(ratio)]);
plotable_mm_small_average_rate_per_user_uplink_ND =  zeros([1 length(ratio)]);
plotable_uhf_small_average_rate_per_user_uplink_ND = zeros([1 length(ratio)]);
plotable_system_coverage_uplink_ND = zeros([1 length(ratio)]);
plotable_tier0_sum_rate_uplink_ND = zeros([1 length(ratio)]);
plotable_tier1_sum_rate_uplink_ND = zeros([1 length(ratio)]);
plotable_tier2_sum_rate_uplink_ND = zeros([1 length(ratio)]);
plotable_tier0_coverage_uplink_ND =  zeros([1 length(ratio)]);
plotable_tier1_coverage_uplink_ND =  zeros([1 length(ratio)]);
plotable_tier2_coverage_uplink_ND =  zeros([1 length(ratio)]);

plotable_system_coverage_ND = zeros([1 length(ratio)]);


plotable_coupled_users_fraction = zeros([1 length(ratio)]);

plotable_decoupled_users_fraction = zeros([1 length(ratio)]);

display('Running Iterations....');

display(strcat('User density:',num2str(user_intensity)));

display(strcat('UHF density:',num2str(UHFM_intensity)));

display(strcat('Downlink Bias UHFM:', num2str(getBiasFactor(0))));
display(strcat('Downlink Bias mmW:', num2str(getBiasFactor(1))));
display(strcat('Downlink Bias UHFS:', num2str(getBiasFactor(2))));
display(strcat('Uplink Bias UHFM:', num2str(getBiasFactorUplink(0))));
display(strcat('Uplink Bias mmW:', num2str(getBiasFactorUplink(1))));
display(strcat('Uplink Bias UHFS:', num2str(getBiasFactorUplink(2))));
display(strcat('Rmin DL:',num2str(R_min_downlink)));
display(strcat('Rmin UL:',num2str(R_min_uplink)));
display('                                                                                                                             ')
tic

for p = 1 : length(ratio)
    display(strcat('Ratio = ', num2str(ratio(p))));
    temp_sumrate_down = 0;      temp_sumrate_up = 0;
    temp_syscov_down  = 0;       temp_syscov_up = 0;
    temp_0_cov_down   = 0;        temp_0_cov_up = 0;
    temp_1_cov_down   = 0;        temp_1_cov_up = 0;
    temp_2_cov_down   = 0;        temp_2_cov_up = 0;
    temp_system_coverage = 0;
    temp_system_coverage_ND = 0;
    
    temp_sumrate_up_ND = 0;
    temp_syscov_up_ND = 0;
    temp_0_cov_up_ND = 0;
    temp_1_cov_up_ND = 0;
    temp_2_cov_up_ND = 0;
    temp_sum_rate_per_tier_uplink_ND = zeros([1 3]);
    temp_uhf_macro_average_rate_per_user_uplink_ND = 0;
    temp_mm_small_average_rate_per_user_uplink_ND  = 0;
    temp_uhf_small_average_rate_per_user_uplink_ND = 0;
    
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
    temp_coupled_users_fraction_ND   = 0;
    temp_decoupled_users_fraction_ND = 0;
    
    for j = 1:Niterations
        
        
        
        %Deploy Basestaions and Users using a PPP
        [nBS_0, nBS_1, nBS_2, BSType, BSLocation] = deployBS(area, UHFM_intensity, ratio(p), mmtoUHFFactor);
        [U, nU] = deployUsers(user_intensity, area);
        
        
        
        
        %Compute angle of mmWave BSs to users
        angle_array = getBsToUserAngle(U, nU, nBS_0, nBS_1, BSLocation);
        
         
        
        
        
        %  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~DOWNLINK BEGIN~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        %%
        %Define downlink variables
        coverage_probability_downlink = 0;
        coverage_matrix_downlink = ones([1 nU]);
        covered_users_per_tier_downlink = zeros([1 3]);
        rate_downlink = zeros([1 nU]);
        
        %Get downlink properties for each Basestation
        [Pt_dB_downlink, B_downlink, f_downlink, subcarriers_downlink] = getBSProperties(BSType);
        
        %Associate Users in downlink
        [nU_BS_0_downlink, nU_BS_1_downlink, nU_BS_2_downlink, Pr_dB_downlink, power_dB_downlink, L_dB_downlink, d, association_matrix_downlink, tier_per_user_downlink, BS_associated_downlink]...
            = downlinkUserAssociation(U, nU, Pt_dB_downlink, nBS_0, nBS_1, nBS_2, BSType, BSLocation);
        
        
        
        %Channels
        [channel, channel_dB] = getChannel(nU, BSType);
        
        %Power control
        %Pt_dB_downlink_new = power_dB_downlink;%
        [Pt_dB_downlink_new] = downlinkPowerControl(U, nU, nBS_0, nBS_1, nBS_2, BSType, BSLocation, association_matrix_downlink, Pt_dB_downlink, Interference_threshold_dB_downlink, channel_dB);
        
        
        
        
        %Get received power
        [Pr_dB_downlink, power_dB_downlink] = getNewReceivedPower(Pt_dB_downlink_new, association_matrix_downlink, U, nU, BSLocation, BSType, nBS_0, nBS_1, nBS_2, channel_dB);
        
        
        %Get bandwidths
        bandwidth_array_downlink = getDownlinkBandwidthPerUser(association_matrix_downlink, BSType, BS_associated_downlink);
        
        %Get downlink SINR
        [SINR_downlink] = downlinkSINR(bandwidth_array_downlink, nU, power_dB_downlink...
            , Pr_dB_downlink, Pt_dB_downlink_new, BS_associated_downlink, BSType, nBS_0, nBS_1, nBS_2, association_matrix_downlink, angle_array, beam_width );
        
        
        %Get downlink Rates
        rate_downlink = bandwidth_array_downlink.* log2(1+(10.^(SINR_downlink/10))).*coverage_matrix_downlink;
        
        %FOR DEBUG
%         BS_associated_downlink;
%         for f = 1 : nU
%         BSTYPE_USER(f) =  BSType(BS_associated_downlink(f));
%         end
        % 
        
        
        %Check Coverage by comparing the value of each user to the set
        %threshold
        for f = 1 : nU
            if rate_downlink(f) < R_min_downlink
                coverage_matrix_downlink(f) = 0;
            else
                covered_users_per_tier_downlink(1+BSType(BS_associated_downlink(f))) = covered_users_per_tier_downlink(1+BSType(BS_associated_downlink(f)))+1;
                temp_sum_rate_per_tier_downlink(1+BSType(BS_associated_downlink(f))) = temp_sum_rate_per_tier_downlink(1+BSType(BS_associated_downlink(f))) + rate_downlink(f);
            end
        end
        
        
        coverage_probability_downlink = sum(coverage_matrix_downlink)/nU;
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
        %%
        % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~UPLINK WITH DECOUPLING ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        
        %Variables
        covered_users_per_tier_uplink = zeros([1 3]);
        SINR_uplink = zeros([1 nU]);
        rate_uplink = zeros([1 nU]);
        coverage_matrix_uplink = ones([1 nU]);
        covered_users_uplink = 0;
        coverage_probability_uplink = 0;
        
        
        P_max_uplink = getUplinkTransmitPower();
        Pt_dB_uplink = P_max_uplink*ones(1, nU);
        
        %Uplink Association
        [BS_associated_uplink, nU_BS_0_uplink, nU_BS_1_uplink, nU_BS_2_uplink, association_matrix_uplink,  Pr_dB_uplink, max_power_db ]...
            = uplinkUserAssociation( Pt_dB_uplink, U, nU, nBS_0, nBS_1, nBS_2, BSType, BSLocation);
        
        %Power Control
        Pt_dB_uplink_new = Pt_dB_uplink;%
        [Pt_dB_uplink_new] = uplinkPowerControl(U, nU, nBS_0, nBS_1, nBS_2, BSType, BSLocation, association_matrix_uplink, Pt_dB_uplink, Interference_threshold_dB_uplink, channel_dB);
        
        %Get received power
        [Pr_dB_uplink_new, power_dB_uplink]...
            = getNewReceivedPower(transpose(Pt_dB_uplink_new), association_matrix_uplink, U, nU, BSLocation, BSType, nBS_0, nBS_1, nBS_2, channel_dB);
        
        %Get Bandwidth per user (It is okay to use downlink function)
        bandwidth_array_uplink = getDownlinkBandwidthPerUser(association_matrix_uplink, BSType, BS_associated_uplink);
        
        %Get interference
        [~, ~, IplusN_db] = getUplinkInterference( BS_associated_uplink,bandwidth_array_uplink,association_matrix_uplink,BSType, Pr_dB_uplink_new, power_dB_uplink, nU);
        
        
        %Get uplink SINR
        SINR_uplink = (power_dB_uplink) - IplusN_db;
        
        %Get uplink Rate
        rate_uplink =  bandwidth_array_uplink.* log2(1+(10.^(SINR_uplink/10))).*coverage_matrix_uplink;
        
        %Get coverage
        for Us = 1:nU
            if rate_uplink(Us) < R_min_uplink
                coverage_matrix_uplink(Us) = 0;
            else
                covered_users_per_tier_uplink(1+BSType(BS_associated_uplink(Us))) = covered_users_per_tier_uplink(1+BSType(BS_associated_uplink(Us))) + 1;
                temp_sum_rate_per_tier_uplink(1+BSType(BS_associated_uplink(Us))) = temp_sum_rate_per_tier_uplink(1+BSType(BS_associated_uplink(Us))) + rate_uplink(Us);
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
        
        %%
         % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~UPLINK WITHOUT DECOUPLING ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        %Variables
        covered_users_per_tier_uplink_ND = zeros([1 3]);
        SINR_uplink_ND = zeros([1 nU]);
        rate_uplink_ND = zeros([1 nU]);
        coverage_matrix_uplink_ND = ones([1 nU]);
        covered_users_uplink_ND = 0;
        coverage_probability_uplink_ND = 0;
        
        
        P_max_uplink_ND = getUplinkTransmitPower();
        Pt_dB_uplink_ND = P_max_uplink_ND*ones(1, nU);
        
        %Uplink Association
        BS_associated_uplink_ND = BS_associated_downlink;
        nU_BS_0_uplink_ND = nU_BS_0_downlink;
        nU_BS_1_uplink_ND = nU_BS_1_downlink;
        nU_BS_2_uplink_ND = nU_BS_2_downlink;
        association_matrix_uplink_ND = association_matrix_downlink;
        
        
        %Power Control
        Pt_dB_uplink_new_ND = Pt_dB_uplink_ND;%
        [Pt_dB_uplink_new_ND] = uplinkPowerControl(U, nU, nBS_0, nBS_1, nBS_2, BSType, BSLocation, association_matrix_uplink_ND, Pt_dB_uplink_ND, Interference_threshold_dB_uplink, channel_dB);
        
        
        
        
        %Get received power
        [Pr_dB_uplink_new_ND, power_dB_uplink_ND]...
            = getNewReceivedPower(transpose(Pt_dB_uplink_new_ND), association_matrix_uplink_ND, U, nU, BSLocation, BSType, nBS_0, nBS_1, nBS_2, channel_dB);
        
        %Get Bandwidth per user (It is okay to use downlink function)
        bandwidth_array_uplink_ND = getDownlinkBandwidthPerUser(association_matrix_uplink_ND, BSType, BS_associated_uplink_ND);
        
        %Get interference
        [~, ~, IplusN_db_ND] = getUplinkInterference( BS_associated_uplink_ND,bandwidth_array_uplink_ND,association_matrix_uplink_ND,BSType, Pr_dB_uplink_new_ND, power_dB_uplink_ND, nU);
        
        
        %Get uplink SINR
        SINR_uplink_ND = (power_dB_uplink_ND) - IplusN_db_ND;
        
        %Get uplink Rate
        rate_uplink_ND =  bandwidth_array_uplink_ND.* log2(1+(10.^(SINR_uplink_ND/10))).*coverage_matrix_uplink_ND;
        
        %Get coverage
        for Us = 1:nU
            if rate_uplink_ND(Us) < R_min_uplink
                coverage_matrix_uplink_ND(Us) = 0;
            else
                covered_users_per_tier_uplink_ND(1+BSType(BS_associated_uplink_ND(Us))) = covered_users_per_tier_uplink_ND(1+BSType(BS_associated_uplink_ND(Us))) + 1;
                temp_sum_rate_per_tier_uplink_ND(1+BSType(BS_associated_uplink_ND(Us))) = temp_sum_rate_per_tier_uplink_ND(1+BSType(BS_associated_uplink_ND(Us))) + rate_uplink_ND(Us);
            end
        end
        coverage_probability_uplink_ND = sum(coverage_matrix_uplink_ND/nU);
        
        
        system_rate_uplink_ND = sum(rate_uplink_ND.*coverage_matrix_uplink_ND);
        
        
        
        if (covered_users_per_tier_uplink_ND(1) ~= 0)
            temp_uhf_macro_average_rate_per_user_uplink_ND = temp_uhf_macro_average_rate_per_user_uplink_ND + (temp_sum_rate_per_tier_uplink_ND(1))./(covered_users_per_tier_uplink_ND(1));
        end
        if (covered_users_per_tier_uplink_ND(2) ~= 0)
            temp_mm_small_average_rate_per_user_uplink_ND = temp_mm_small_average_rate_per_user_uplink_ND + (temp_sum_rate_per_tier_uplink_ND(2))./(covered_users_per_tier_uplink_ND(2));
        end
        if (covered_users_per_tier_uplink_ND(3) ~= 0)
            temp_uhf_small_average_rate_per_user_uplink_ND = temp_uhf_small_average_rate_per_user_uplink_ND + (temp_sum_rate_per_tier_uplink_ND(3))./(covered_users_per_tier_uplink_ND(3));
        end
        
        temp_sumrate_up_ND = temp_sumrate_up_ND + system_rate_uplink_ND;
        temp_syscov_up_ND = temp_syscov_up_ND + coverage_probability_uplink_ND ;
        temp_0_cov_up_ND = temp_0_cov_up_ND + (covered_users_per_tier_uplink_ND(1)/nU) ;
        temp_1_cov_up_ND = temp_1_cov_up_ND + (covered_users_per_tier_uplink_ND(2)/nU) ;
        temp_2_cov_up_ND = temp_2_cov_up_ND + (covered_users_per_tier_uplink_ND(3)/nU) ;
        temp_system_coverage_ND = temp_system_coverage_ND + (sum(coverage_matrix_downlink.*coverage_matrix_uplink_ND)/nU);
        
    end
    clc
    display(strcat('Percent Completion:  ', strcat(num2str(ceil(100*p/length(ratio)))),'%'));

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
    
    
    
    
    
    
    %----------
    plotable_system_coverage_ND(p) = temp_system_coverage_ND/Niterations;
    plotable_uhf_macro_average_rate_per_user_uplink_ND(p) = temp_uhf_macro_average_rate_per_user_uplink_ND/Niterations;
    plotable_mm_small_average_rate_per_user_uplink_ND(p) = temp_mm_small_average_rate_per_user_uplink_ND/Niterations;
    plotable_uhf_small_average_rate_per_user_uplink_ND(p) = temp_uhf_small_average_rate_per_user_uplink_ND/Niterations;
    plotable_sum_rate_uplink_ND(p) = temp_sumrate_up_ND/Niterations;
    plotable_system_coverage_uplink_ND(p) = temp_syscov_up_ND/Niterations;
    plotable_tier0_sum_rate_uplink_ND(p) = temp_sum_rate_per_tier_uplink_ND(1)/Niterations;
    plotable_tier1_sum_rate_uplink_ND(p) = temp_sum_rate_per_tier_uplink_ND(2)/Niterations;
    plotable_tier2_sum_rate_uplink_ND(p) = temp_sum_rate_per_tier_uplink_ND(3)/Niterations;
    
    plotable_tier0_coverage_uplink_ND(p) = temp_0_cov_up_ND/Niterations;
    plotable_tier1_coverage_uplink_ND(p) = temp_1_cov_up_ND/Niterations;
    plotable_tier2_coverage_uplink_ND(p) = temp_2_cov_up_ND/Niterations;
    %----------
    
    
    
    
    
    plotable_coupled_users_fraction(p) = temp_coupled_users_fraction/Niterations;
    plotable_decoupled_users_fraction(p) = temp_decoupled_users_fraction/Niterations;
end

display(strcat('Saved to ', general_filename))
    
if PC
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
hold on
plot(ratio, plotable_tier0_sum_rate_downlink, 'g');
plot(ratio, plotable_tier1_sum_rate_downlink, 'k');
plot(ratio, plotable_tier2_sum_rate_downlink, 'm');
legend('Tier 0 Sumrate', 'Tier 1 sumrate', 'Tier 2 sumrate');



figure;
hold on
plot(ratio, plotable_tier0_coverage_uplink, 'r--*');
title('Coverage Probability - Uplink ');
xlabel('l_s / l_u_h_f');
ylabel('Coverage Probability');
plot(ratio, plotable_tier1_coverage_uplink, 'g--^');
plot(ratio, plotable_tier2_coverage_uplink, 'k--<');
plot(ratio, plotable_system_coverage_uplink, 'm-->');
plot(ratio, plotable_tier0_coverage_uplink_ND, 'r:s');
plot(ratio, plotable_tier1_coverage_uplink_ND, 'g:d');
plot(ratio, plotable_tier2_coverage_uplink_ND, 'k:h');
plot(ratio, plotable_system_coverage_uplink_ND, 'm:p');
legend('UHF Macro', 'mmWave Small', 'UHF Small', 'Overall', 'UHF Macro ND', 'mmWave Small ND', 'UHF Small ND', 'Overall ND');

figure;
plot(ratio, plotable_system_coverage, 'm--s');hold on
plot(ratio, plotable_system_coverage_ND, 'm--d');
title('System Coverage');
xlabel('l_s / l_u_h_f');
legend('Decoupled','Non decoupled')
figure;
hold on
plot(ratio, plotable_sum_rate_uplink, 'r-s'); hold on
plot(ratio, plotable_sum_rate_uplink_ND, 'g-d');
legend('Decoupled', 'Not decoupled')
title('Rate - Uplink');
xlabel('l_s / l_u_h_f');
ylabel('Rate(bps)');
figure;
hold on
plot(ratio, plotable_tier0_sum_rate_uplink, 'g');
hold on
plot(ratio, plotable_tier1_sum_rate_uplink, 'k');
plot(ratio, plotable_tier2_sum_rate_uplink, 'm');
plot(ratio, plotable_tier0_sum_rate_uplink_ND, 'g--d');
plot(ratio, plotable_tier1_sum_rate_uplink_ND, 'k--*');
plot(ratio, plotable_tier2_sum_rate_uplink_ND, 'm--s');
title('Rate per user - Uplink');
xlabel('l_s / l_u_h_f');
ylabel('Rate(bps)');
legend('Tier 0 sum rate', 'Tier 1 sum rate', 'Tier 2 sum rate', 'Tier 0 sum rate ND', 'Tier 1 sum rate ND', 'Tier 2 sum rate ND');


figure;
plot(ratio, plotable_coupled_users_fraction, ratio, plotable_decoupled_users_fraction);
title('Decoupling');
xlabel('l_s / l_u_h_f');
ylabel('Fraction of users');
legend('Coupled', 'Decoupled');
display('Done');
end

general_filename = strcat(general_filename,seperator);

save(strcat(general_filename, 'ratio.mat'), 'ratio');
% saving downlink variables
save(strcat(general_filename, 'plotable_uhf_macro_average_rate_per_user_downlink.mat'), 'plotable_uhf_macro_average_rate_per_user_downlink');
save(strcat(general_filename, 'plotable_mm_small_average_rate_per_user_downlink.mat'), 'plotable_mm_small_average_rate_per_user_downlink');
save(strcat(general_filename, 'plotable_uhf_small_average_rate_per_user_downlink.mat'), 'plotable_uhf_small_average_rate_per_user_downlink');
save(strcat(general_filename, 'plotable_sum_rate_downlink.mat'), 'plotable_sum_rate_downlink');
save(strcat(general_filename, 'plotable_tier0_sum_rate_downlink.mat'), 'plotable_tier0_sum_rate_downlink');
save(strcat(general_filename, 'plotable_tier1_sum_rate_downlink.mat'), 'plotable_tier1_sum_rate_downlink');
save(strcat(general_filename, 'plotable_tier2_sum_rate_downlink.mat'), 'plotable_tier2_sum_rate_downlink');
save(strcat(general_filename, 'plotable_system_coverage_downlink.mat'), 'plotable_system_coverage_downlink');
save(strcat(general_filename, 'plotable_tier0_coverage_downlink.mat'), 'plotable_tier0_coverage_downlink');
save(strcat(general_filename, 'plotable_tier1_coverage_downlink.mat'), 'plotable_tier1_coverage_downlink');
save(strcat(general_filename, 'plotable_tier2_coverage_downlink.mat'), 'plotable_tier2_coverage_downlink');

save(strcat(general_filename, 'plotable_system_coverage.mat'), 'plotable_system_coverage');
%saving uplink variables

save(strcat(general_filename, 'plotable_uhf_macro_average_rate_per_user_uplink.mat'), 'plotable_uhf_macro_average_rate_per_user_uplink');
save(strcat(general_filename, 'plotable_mm_small_average_rate_per_user_uplink.mat'), 'plotable_mm_small_average_rate_per_user_uplink');
save(strcat(general_filename, 'plotable_uhf_small_average_rate_per_user_uplink.mat'), 'plotable_uhf_small_average_rate_per_user_uplink');
save(strcat(general_filename, 'plotable_sum_rate_uplink.mat'), 'plotable_sum_rate_uplink');
save(strcat(general_filename, 'plotable_tier0_sum_rate_uplink.mat'), 'plotable_tier0_sum_rate_uplink');
save(strcat(general_filename, 'plotable_tier1_sum_rate_uplink.mat'), 'plotable_tier1_sum_rate_uplink');
save(strcat(general_filename, 'plotable_tier2_sum_rate_uplink.mat'), 'plotable_tier2_sum_rate_uplink');
save(strcat(general_filename, 'plotable_system_coverage_uplink.mat'), 'plotable_system_coverage_uplink');
save(strcat(general_filename, 'plotable_tier0_coverage_uplink.mat'), 'plotable_tier0_coverage_uplink');
save(strcat(general_filename, 'plotable_tier1_coverage_uplink.mat'), 'plotable_tier1_coverage_uplink');
save(strcat(general_filename, 'plotable_tier2_coverage_uplink.mat'), 'plotable_tier2_coverage_uplink');

%---------
save(strcat(general_filename, 'plotable_uhf_macro_average_rate_per_user_uplink_ND.mat'), 'plotable_uhf_macro_average_rate_per_user_uplink_ND');
save(strcat(general_filename, 'plotable_mm_small_average_rate_per_user_uplink_ND.mat'), 'plotable_mm_small_average_rate_per_user_uplink_ND');
save(strcat(general_filename, 'plotable_uhf_small_average_rate_per_user_uplink_ND.mat'), 'plotable_uhf_small_average_rate_per_user_uplink_ND');
save(strcat(general_filename, 'plotable_sum_rate_uplink_ND.mat'), 'plotable_sum_rate_uplink_ND');
save(strcat(general_filename, 'plotable_tier0_sum_rate_uplink_ND.mat'), 'plotable_tier0_sum_rate_uplink_ND');
save(strcat(general_filename, 'plotable_tier1_sum_rate_uplink_ND.mat'), 'plotable_tier1_sum_rate_uplink_ND');
save(strcat(general_filename, 'plotable_tier2_sum_rate_uplink_ND.mat'), 'plotable_tier2_sum_rate_uplink_ND');
save(strcat(general_filename, 'plotable_system_coverage_uplink_ND.mat'), 'plotable_system_coverage_uplink_ND');
save(strcat(general_filename, 'plotable_tier0_coverage_uplink_ND.mat'), 'plotable_tier0_coverage_uplink_ND');
save(strcat(general_filename, 'plotable_tier1_coverage_uplink_ND.mat'), 'plotable_tier1_coverage_uplink_ND');
save(strcat(general_filename, 'plotable_tier2_coverage_uplink_ND.mat'), 'plotable_tier2_coverage_uplink_ND');
save(strcat(general_filename, 'plotable_system_coverage_ND.mat'), 'plotable_system_coverage_ND');
save(strcat(general_filename, 'plotable_system_coverage.mat'), 'plotable_system_coverage');

%---------------

%decoupling variables
save(strcat(general_filename, 'plotable_coupled_users_fraction.mat'), 'plotable_coupled_users_fraction');
save(strcat(general_filename, 'plotable_decoupled_users_fraction.mat'), 'plotable_decoupled_users_fraction');

save(strcat(general_filename, 'Niterations.mat'), 'Niterations');





display(strcat('Task finished in ', strcat(num2str(toc/60), ' minutes.')))