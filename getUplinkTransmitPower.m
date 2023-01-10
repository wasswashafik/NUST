function [ Pt_dB ] = getUplinkTransmitPower()
%UNTITLED10 Summary of this function goes here
%   Detailed explanation goes here


%CHANGE POWER TO 20dBm
POWER = 24; %power in dBm
Pt_dB = POWER - 30;

end
