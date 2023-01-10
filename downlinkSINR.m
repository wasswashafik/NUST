function [SINR_downlink] = downlinkSINR(bandwidth_array_downlink, nU, power_dB_downlink, Pr_dB_downlink, Pt_dB_downlink, BS_associated_downlink, BSType, nBS_0, nBS_1, nBS_2, association_matrix_downlink, angle_array, beam_width , SINR_Threshold)

        SINR_downlink = zeros([1, nU]);
       

         % AFTER ASSOCIATION, GET DOWNLINK INTERFERENCE AT EACH USER
         [interference_watts_downlink, interference_dB_downlink] = ...
             getDownlinkInterference(BS_associated_downlink, BSType, nBS_0, nBS_1, nBS_2, association_matrix_downlink, Pr_dB_downlink, angle_array, beam_width);
%         
         %GET BANDWIDTH THAT EACH USER HAS ( BW OF BS / NUMBER OF USER ON
         %THAT BS)
         IplusN_downlink = getDownlinkIplusN(interference_watts_downlink, bandwidth_array_downlink);
         SINR_downlink = (power_dB_downlink) - IplusN_downlink;
 
             

end