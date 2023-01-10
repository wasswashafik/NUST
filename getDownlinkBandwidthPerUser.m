function [bandwidth_array] = getDownlinkBandwidthPerUser(association_matrix, BSType, BS_associated)

% BS_associated -> returns number of BS associated to each user, eg
% BS_associated(1) = 2 means that user 1 is connected to BS 2

%BSType(2) = 0 means that BS 2 is tier 0. There are 3 tiers, 0, 1, and 2

%association_matrix = binary array showing which user is connected to which
%BS

nU = length(association_matrix(:, 1));
nBS = length(association_matrix(1, :));
bandwidth_array = zeros([1 nU]);
user_per_BS = zeros([1 nBS]);
bandwidth_per_BS = zeros([1 length(association_matrix(1, :))]);

for BS = 1:length(association_matrix(1, :)) % for each BS
    
    user_per_BS(BS) = sum(association_matrix(:, BS));
    bandwidth_per_BS(BS) = getBandwidth(BSType(BS))/(2.*user_per_BS(BS));       %half bandwidth allocated to downlink
    
end
 
for user = 1:nU
    bandwidth_array(user) = bandwidth_per_BS(BS_associated(user));
end

end