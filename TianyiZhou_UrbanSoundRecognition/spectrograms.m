% Zach Neveu | Christian Grenier | Tianyi Zhou
%
% Compute spectrograms for all files in an audio dataset and return the resulting data.
% mostly from 
% https://www.mathworks.com/help/deeplearning/examples/deep-learning-speech-recognition.html?s_tid=mwa_osa_a 
% with some changes by Zach (single channel)

function X = spectrograms(ads,segmentDuration,frameDuration,hopDuration,numBands)

numHops = ceil((segmentDuration - frameDuration)/hopDuration);
numFiles = length(ads.Files);
X = zeros([numBands,numHops+1,1,numFiles],'single');

for i = 1:numFiles

	% Read an item from the dataset that was passed
	[x,info] = read(ads);

	% get function params from data given
	fs = info.SampleRate;
	frameLength = round(frameDuration*fs);
	hopLength = round(hopDuration*fs);

	try
		spec = melSpectrogram(x,fs, ...
		'WindowLength',frameLength, ...
		'OverlapLength',frameLength - hopLength, ...
		'NumBands',numBands, ...
		'FrequencyRange',[50,20000]);

		% For samples less than 4 seconds, place the spectrogram in the
		% middle of the space
		width = size(spec,2);
		left = floor((numHops-width)/2)+1;
		ind = left+1:left+width;

		% Take only first channel
		X(:,ind,1,i) = spec(:,:,1);

	catch e
		disp("Failed for file "+i)
		disp(info)
		disp(e)
		disp(size(x))
		continue
	end


	if mod(i,100) == 0
		disp("Processed " + i + " files out of " + numFiles)
	end

end

end
