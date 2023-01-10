function [LoS] = isLOS(link_distance)
% This funtion returns binary values, 1 indicating LOS and 0 indicating
% NLOS link for mmWave links. Blockages are modelled using a Boolean model
% of rectangles based on the random shape theory as proposed in the
% following literature:
% T. Bai, R. Vaze, R. W. Heath, “Analysis of blockage effects on urban
% cellular networks,” IEEE Trans. Wireless Commun., vol. 13, no. 9, pp.
% 5070–5083, Sept. 2014.

beta = 0.0014;          % Beta parameter for NUST Campus beta = 0.0014;

R = link_distance;
p_LoS = exp(-beta*R);   % Generates a PDF for the link being LoS

% Pseudorandomly chooses if the link is LoS or NLoS based on p_LoS
if rand(1) <= p_LoS
   LoS = 1;
else
   LoS = 0;
end

end