function im_scrambled = phase_scramble(im,showimg)
    if nargin < 2
        showimg = 0;
    end
    % Convert the image to grayscale
    im_gray = rgb2gray(im);
    % Compute the Fourier transform of the image
    im_fft = fft2(double(im_gray));
    % Extract the magnitude and phase information
    im_mag = abs(im_fft);
    im_phase = angle(im_fft);
    % Randomize the phase information
    [m,n] = size(im_phase);
%     rand_phase = -pi + (pi+pi)*rand(size(im_phase)); % Between -pi and pi
%     im_phase_scrambled = im_phase.*exp(1i*rand_phase);
    im_phase_scrambled = reshape(im_phase(randperm(numel(im_phase))),[m n]);
    % Combine the randomized phase information with the original magnitude information
    im_scrambled_fft = im_mag.*exp(1i*im_phase_scrambled);
    % Compute the inverse Fourier transform to obtain the phase-scrambled image
    im_scrambled = ifft2(im_scrambled_fft);
    % Convert the phase-scrambled image back to uint8 data type
    im_scrambled = uint8(abs(im_scrambled));
    % Display the original and phase-scrambled images side by side
    if showimg
        subplot(1,2,1), imshow(im_gray), title('Original');
        subplot(1,2,2), imshow(im_scrambled), title('Phase-scrambled');
    end
end