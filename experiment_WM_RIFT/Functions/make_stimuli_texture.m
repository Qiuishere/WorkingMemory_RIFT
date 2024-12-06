   
function [imgTxts, maskTxts, prm] = make_stimuli_texture(prm)

%% targets
    imgTxts = zeros(length(prm.fac.figures), length(prm.fac.adjustRange), length(prm.fac.views));
    for i = 1:length(prm.fac.figures)
        for j = prm.fac.adjustRange(1): prm.fac.adjustRange(end)
            for k = 1:length(prm.fac.views)
                theimg = strcat('stimuli/Front-Upper-shoulder-middle/', prm.fac.figures{i}, '_Front_', num2str(j), prm.fac.views{k}, '.png');
                [X1, ~, alpha1]  = imread(theimg);
                prm.img.scale    = prm.img.WPix/size(X1,2);
                
                X1               = imresize(X1,prm.img.scale);
                alpha1           = imresize(alpha1, prm.img.scale);
                
                X1(:,:,4) = alpha1;
                imgTxts(i,j,k) =   Screen( 'MakeTexture', prm.w.Number, X1);
            end
        end
    end
    
    prm.img.W      = size(X1,2);
    prm.img.H      = size(X1,1);
    prm.img.offPix = 139 * prm.img.scale;% distance from the shoulder to the center of the image. vertical shift to center the shoulder at fixation

    prm.img.rect   = [0, 0, prm.img.W, prm.img.W];
    prm.img.presentedSize = [ min(prm.img.W, prm.w.Center(1)), min(prm.img.H, prm.w.Center(2))];
    prm.img.sourceRect = CenterRectOnPoint([0, 0, prm.img.presentedSize], prm.img.W/2, prm.img.H/2 - prm.img.offPix);
    prm.img.sourceRect(prm.img.sourceRect<0) = 0;
        
    prm.bar.armLength = 216 * prm.img.scale;
    prm.bar.color = 0.6*255; % a light grey
    
%% mask
prm.exp.Nmask = 6;
    for theimg = 1:prm.exp.Nmask
        maskFile = strcat('stimuli/masks/m', num2str(theimg), '.jpg');
        mask = imread(maskFile);
        mask = imresize(mask, prm.img.WPix/size(mask,2));
        maskTxts(theimg)  =   Screen( 'MakeTexture', prm.w.Number, mask);
    end