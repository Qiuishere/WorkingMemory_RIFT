
function [blackTxts, whiteTxts] = make_tag_texture(prm)

imgSize = RectWidth(prm.img.pos) ;

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
sigma = 2; % Adjust sigma for desired smoothness
%img = imgaussfilt(img, sigma);
%img(Y,X) = 255;
alphaLayer = img;

figure;
imshow(img/255)


for i = 1:prm.tag.nStrip
    theimg = imrotate(img, prm.tag.angle(i), 'bilinear', 'crop');
    thealpha = imrotate(alphaLayer, prm.tag.angle(i), 'bilinear', 'crop'); % 'bilinear' interpolation, 'crop' option
    blackbar = cat(3, 255-theimg, thealpha);
    whitebar = cat(3, theimg, thealpha);
    
    
    blackTxts(i,1) = Screen('MakeTexture', prm.w.Number,  blackbar);
    blackTxts(i,2) = Screen('MakeTexture', prm.w.Number,  fliplr(blackbar));
    whiteTxts(i,1) = Screen('MakeTexture', prm.w.Number,  whitebar);
    whiteTxts(i,2) = Screen('MakeTexture', prm.w.Number,  fliplr(whitebar));
end



