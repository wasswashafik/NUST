function [bias_factor_db] = getUplinkAntennaGain(tier, LOS_NLOS)
%This function returns the biasing factor for the concerned tier

if tier == 0
    bias_factor_db = 0;  %No biasing for UHF macrocell
elseif tier == 1
    if LOS_NLOS == 1
        bias_factor_db = 18; %mmWave small cell biased by 20dB
    else
        bias_factor_db = 0; %mmWave small cell biased by 20dB
    end
    
elseif tier == 2
    bias_factor_db=0;    %uhf small cell bias factor
end
end