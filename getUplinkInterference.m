function [uplink_interference, uplink_interference_dB, IplusN_db] = getUplinkInterference( BS_associated,bandwidth_array,association_matrix_uplink,BSType, Pr_dB, max_power_db, nU)
nBS = length(association_matrix_uplink(1, :));
uplink_interference = zeros([1 nU ]);
%For each User, 

% -> For all Base Stations other than the one concerned
% -> Choose a random user associated with it
% -> add that users power as an interference to the concerned BS

power_received_watts = 10.^(Pr_dB./10); %uplink received power array
noise = zeros([1 nU]);


for concerned_user = 1:nU
    
   [A,noise(concerned_user)] = getNoise(bandwidth_array(concerned_user));
    
     for interfering_user = 1:nU
         if concerned_user == interfering_user
             continue
         end
         if BSType( BS_associated(concerned_user)) ~= BSType(BS_associated(interfering_user))
             continue
         end
    
     uplink_interference(concerned_user) = uplink_interference(concerned_user)+power_received_watts(interfering_user, BS_associated(concerned_user));
     [A,noise(concerned_user)] = getNoise(bandwidth_array(concerned_user));
    
    end
end

IplusN_db = 10.*log10( uplink_interference + noise );



uplink_interference_dB = 10.*log10(uplink_interference );




end