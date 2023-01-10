function [noise_dB, noise_watts] = getNoise(BW)
%funtion to return noise.
%values obtianed from "Self-adaptive Power Control Mechanism in
%Device-to-Device Enabled Multi-tier HetNets; An
%Optimization Approach"


N0 = -174-30; %noise in dB/Hz
N0 = 10.^(N0./10); %noise converted to Watts/Hz
noise_watts = BW*N0;
noise_dB = 10.*log10(noise_watts);


end