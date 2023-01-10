function [IplusN] = getDownlinkIplusN(interference_watts, bandwidth_array)
nU = length(bandwidth_array);  %This length is same as number of users as bandwidth array gives the bandwidth for each user
IplusN = zeros([1 nU]);
for user = 1:nU
     [noise_dB, noise_watts] = getNoise(bandwidth_array(user));
    IplusN(user) =  interference_watts(user) + noise_watts;
    
end

IplusN = 10.*log10(IplusN);

end