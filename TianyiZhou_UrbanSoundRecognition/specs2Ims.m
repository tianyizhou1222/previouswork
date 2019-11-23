% Zach Neveu | Tianyi Zhou | Christian Grenier
%
% Compute color images from magnitude spectrograms
%
% When using pre-trained networks, the input sizes are fixed, as well as the number of channels.
% This function allows for converting spectrograms generated with spectrograms.m into valid
% input images.


function Y = specs2Ims(X, sz)
	XIms = zeros(sz(1), sz(2), 3, size(X, 4));
	for i=1:size(X, 4)
		close all
		figure('visible', 'off')
		pcolor(X(:,:,i));
		shading flat
		Frame = getframe;
		[imi, ~] = frame2im(Frame);
		XIms(:,:,:,i) = imresize(imi, [sz(1) sz(2)]);
		if mod(i, 50) == 0
			disp("Processed "+num2str(i))
		end
	end
	Y = XIms;
end

