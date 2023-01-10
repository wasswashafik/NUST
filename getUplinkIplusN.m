function [IplusN_uplink] = getUplinkIplusN(interference_watts_uplink, bandwidth_array_uplink)

nU = length(bandwidth_array_uplink);  %This length is same as number of users as bandwidth array gives the bandwidth for each user
IplusN_uplink = zeros([1 nU]);
for user = 1:nU
     [noise_dB, noise_watts] = getNoise(bandwidth_array_uplink(user));
    IplusN_uplink(user) =  interference_watts_uplink(user) + noise_watts;
    
end

IplusN_uplink = 10.*log10(IplusN_uplink);

end