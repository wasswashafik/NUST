function [fade_channel, fade_channel_dB ] = getChannel(nU, BSType)

fade_channel = zeros(nU, length(BSType));

for BS = 1:length(BSType)
    fade_channel(:, BS) = fading(nU, BSType(BS));
end

fade_channel_dB = 10.*log10(fade_channel);

end
