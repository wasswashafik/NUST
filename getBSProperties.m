function [ Pt_dB, B, f, subcarriers ] = getBSProperties( tier )
Pt_dB = getTransmitPower(tier);
B = getBandwidth(tier);
f = getFrequency(tier);
subcarriers = getSubCarriers(tier);
end

