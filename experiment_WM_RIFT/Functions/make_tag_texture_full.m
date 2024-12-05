
function [blackTxts] = make_tag_texture(prm)

imgSize = prm.img.WPix ;

armLength = round(prm.bar.armLength);

img = zeros(imgSize, imgSize);

if mod(imgSize, 2)==0
    center = imgSize/2;
else
    center = ceil(imgSize/2);
end
% Create meshgrid of pixel coordinates relative to center
[xGrid, yGrid] = meshgrid(1:imgSize, 1:imgSize);
x = xGrid - center;
y = center - yGrid; % Flip y to make it in the "upright" direction

% Convert Cartesian coordinates to polar coordinates
[theta, radius] = cart2pol(x, y);
theta = rad2deg(theta); % Convert theta to degrees

% Define the fan angle and radius range
fanAngle = prm.tag.stripWidth; % Total angle of fan in degrees
halfFanAngle = fanAngle / 2; % Half of the angle for symmetric bounds
angleCenter = 90; % Upright direction (90 degrees)
minAngle = angleCenter - halfFanAngle;
maxAngle = angleCenter + halfFanAngle;

% Define radius range for the fan
minRadius = armLength - round(0.5 * armLength);
maxRadius = armLength;

% Create a mask for the fan area
fanMask = (radius >= minRadius & radius <= maxRadius) & ...
    (theta >= minAngle & theta <= maxAngle);

% Set pixels within the fan mask to 255 (white)
img(fanMask) = 255;


% Apply Gaussian smoothing
%sigma = 2; % Adjust sigma for desired smoothness
%img = imgaussfilt(img, sigma);
%img(Y,X) = 255;
alphaLayer = img;

figure

blackTxts = zeros(2, size(prm.tag.tag_sigs, 2));
for thefr = 1: 30%size(prm.tag.tag_sigs, 2)
    for i = 1:prm.tag.nStrip
        
        thelumi  = prm.tag.tag_sigs(i, thefr);
        theimg   = img *thelumi;
        fans(:,:,i)   = imrotate(theimg, prm.tag.angle(i), 'bilinear', 'crop');
        % thealpha = imrotate(alphaLayer, prm.tag.angle(i), 'bilinear', 'crop'); % 'bilinear' interpolation, 'crop' option
        
    end
    
    allfans = sum(fans,3);
    
    subplot(3,10,thefr)
    imshow(allfans/255)
    alpha = zeros(imgSize, imgSize);
    alpha(allfans~=0) = 255;
    texture = repmat(allfans,[1,1,3]);
    texture(:,:,4) = alpha;
    
    
    blackTxts(thefr,1) = Screen('MakeTexture', prm.w.Number,  texture);
    %blackTxts(thefr,2) = Screen('MakeTexture', prm.w.Number,  fliplr(texture));
    
end



