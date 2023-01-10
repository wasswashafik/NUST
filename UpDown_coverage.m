clc;
clear all
close all
% -------------------------------------------------------------------------

area = 1000*1000;                  % the total area in squared meters

user_intensity = 1e-4;          % the intensity of users in m^2
UHFM_intensity = 5e-6;              % the BS intensity in m^2
ratio = 0:20;
mm_intensity = ratio*UHFM_intensity;
%deploy users in the area

Niterations = 10;
warning off;
threshold_dB = -999;
uplink_threshold = -999;
temp_cov_prob_down = 0;
temp_cov_prob_mm_down = 0;
coverage_probability_down = zeros(1, length(ratio));
coverage_probability_mmWave_down = zeros(1, length(ratio));

temp_cov_prob_up = 0;
temp_cov_prob_mm_up = 0;
temp_decoupled_user = 0;
temp_coupled_user = 0;
temp_downlink_macro = 0;
temp_uplink_macro = 0;
temp_downlink_small = 0;
temp_uplink_small = 0;


coverage_probability_up = zeros(1, length(ratio));
coverage_probability_mmWave_up = zeros(1, length(ratio));

downlink_macro = zeros(1 , length(ratio));
uplink_macro = zeros(1 , length(ratio));
downlink_small = zeros(1 , length(ratio));
uplink_small = zeros(1 , length(ratio));
decoupled_user = zeros(1 , length(ratio));
coupled_user = zeros(1 , length(ratio));

for p = 1 : length(ratio)
   

    for j = 1:Niterations
        [nBS_0, nBS_1, nBS_2, BSType, BSLocation] = deployBS(area, UHFM_intensity, ratio(p), 0.1);
        [U, nU] = deployUsers(user_intensity, area);
        covered_users = 0;
        % DOWNLINK
        downlink_user_association;
        temp_cov_prob_down = temp_cov_prob_down + (covered_users/nU);
        temp_cov_prob_mm_down = temp_cov_prob_mm_down + ((nU_BS_1+nU_BS_2)/covered_users);
        % DOWNLINK
        uplink_covered_users = 0;
        % UPLINK
        
        uplink_user_association;
        temp_cov_prob_up = temp_cov_prob_up + (uplink_covered_users/nU);
        temp_cov_prob_mm_up = temp_cov_prob_mm_up + ((uplink_nU_BS_1+uplink_nU_BS_2)/uplink_covered_users);
        temp_coupled_user =  temp_coupled_user +  (sum(sum(association_matrix_uplink.*association_matrix)))./nU;
        temp_decoupled_user =  temp_decoupled_user + ((length(Ux) - ((sum(sum(association_matrix_uplink.*association_matrix))))))./nU;
        % UPLINK
        temp_downlink_macro = temp_downlink_macro + ((sum(sum(association_matrix(:, 1:nBS_0)))))./nU;
        temp_uplink_macro = temp_uplink_macro+ sum(sum(association_matrix_uplink(:, 1:nBS_0)))./nU;
        temp_downlink_small = temp_downlink_small + sum(sum(association_matrix(:, nBS_0+1:(nBS_0+nBS_1+nBS_2))))./nU;
        temp_uplink_small = temp_uplink_small+ sum(sum(association_matrix_uplink(:, nBS_0+1:(nBS_0+nBS_1+nBS_2))))./nU;
        %Values to analyze

        clc;
        display(strcat('Percent Completion:  ', strcat(num2str(ceil(100*p/length(ratio)))),'%'));
        display(strcat('Percent Completion:  ', strcat(num2str(ceil(100*j/Niterations))),'%'));

    end
    % DOWNLINK
    coverage_probability_down(p) = (temp_cov_prob_down/Niterations);
    coverage_probability_mmWave_down(p) = (temp_cov_prob_mm_down/Niterations);
    coupled_user (p) = temp_coupled_user/Niterations;
    decoupled_user (p) = temp_decoupled_user/Niterations;
    downlink_macro(p) = temp_downlink_macro/Niterations;
    uplink_macro(p) = temp_uplink_macro/Niterations;
    downlink_small(p) = temp_downlink_small/Niterations;
    uplink_small(p) = temp_uplink_small/Niterations;
    
    temp_cov_prob_down=0;
    temp_cov_prob_mm_down=0;
    temp_coupled_user = 0;
    temp_decoupled_user = 0;
    % DOWNLINK
    % UPLINK
    coverage_probability_up(p) = (temp_cov_prob_up/Niterations);
    coverage_probability_mmWave_up(p) = (temp_cov_prob_mm_up/Niterations);
    temp_cov_prob_up=0;
    temp_cov_prob_mm_up=0;
    temp_downlink_macro = 0;
    temp_uplink_macro = 0;
    temp_downlink_small = 0;
    temp_uplink_small = 0;
    % UPLINK
   
end
figure;
plot(ratio, coverage_probability_mmWave_down, 'r');
title('Fraction of users attended to small cells ');
xlabel('l_s / l_u_h_f');
ylabel('Fraction of users of small cell ');
hold on
plot(ratio, coverage_probability_mmWave_up, 'g');
legend('DOWNLINK', 'UPLINK');

figure;
plot(ratio, decoupled_user, 'r');
title('Fraction of decoupled and coupled users ');
xlabel('l_s / l_u_h_f');
ylabel('Fraction of users');
hold on
plot(ratio, coupled_user, 'g');
legend('Decoupled', 'Coupled');



figure;
plot(ratio, downlink_macro, 'r');
title('Fraction of user tagged to Macro/Small cell in uplink/downlink ');
xlabel('l_s / l_u_h_f');
ylabel('Fraction of users');
hold on;
plot(ratio, downlink_small, 'g');
hold on;
plot(ratio, uplink_macro, 'k');
hold on;
plot(ratio, uplink_small, 'm');




legend('Macro-downlink', 'Small-downlink', 'Macro-uplink', 'Small-uplink');


% figure;
% plot(ratio, coverage_probability);
% title('Fraction of users attended');
% xlabel('l_s / l_u_h_f');
% ylabel('Fraction of users');





display('Done');
