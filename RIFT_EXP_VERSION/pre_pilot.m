%% RIFT pre-pilot ===============================================
% Songyun at DCC
%                                                                11.01.2024
%% ========================================================================
Screen('Preference', 'SkipSyncTests', 1);
sca;
addpath('/Users/songyunbai/Documents/Donders/Experiment/RIFT test/Matlab/')
prm = [];
% general exp prm
prm.exp.name           = 'RIFT Pre-pilot';
prm.exp.path           = pwd;
prm.exp.datetime       = datetime;

prm.screen.refresh     = 60*12;
%% monitor 
is_live                = 0;
ptb                    = InitPsychtoolbox(is_live);


%% stimuli prm
prm.stim.path          = "/Users/songyunbai/Documents/Donders/Experiment/RIFT test/dataset/OrigBigAnimals/";
prm.stim.screen_w      = 48;%53;%38;     % cm (38 2nd screen at desk; 53 in cubicle or home desk; 48 in meg)
prm.stim.view_dist     = 85;%61;%70;    % cm (70 2nd screen at desk; 61 in cubicle or home desk (approx); 85 in meg)
prm.stim.pixel_size    = prm.stim.screen_w / ptb.win_w * 2; % note: *2 factor because of ProPIXX tagging
prm.stim.pix2deg       = @(px) (360./pi .* atan(px.*prm.stim.pixel_size./(2.*prm.stim.view_dist)));
prm.stim.deg2pix       = @(deg) (2.*prm.stim.view_dist.*tan(pi./360.*deg)./prm.stim.pixel_size);
prm.stim.tag_frex      = [5,10,15,20];
prm.stim.duration      = 10; % in second
prm.stim.dt            = 1 / prm.screen.refresh; % in second
prm.stim.timax         = prm.stim.dt:prm.stim.dt:prm.stim.duration;
prm.stim.tag_sigs      = cos(2*pi*prm.stim.tag_frex(:)*prm.stim.timax) / 2 + 0.5;
prm.stim.size          = [100,100];
prm.stim.subcenters    = [-ptb.win_w/4, -ptb.win_h/4; ptb.win_w/4, -ptb.win_h/4;...
                          -ptb.win_w/4, ptb.win_h/4; ptb.win_w/4, ptb.win_h/4];
prm.stim.distance      = [-ptb.win_w/16, -ptb.win_h/16; ptb.win_w/16, -ptb.win_h/16;...
                          -ptb.win_w/16, ptb.win_h/16; ptb.win_w/16, ptb.win_h/16];
prm.stim.disControl    = [-1,0,1,1];               
                      
                      
cd(prm.stim.path)
ss              = indir;
n               = numel(ss);
selected_img    = randsample(n,4) - 1;

Texture_img = zeros(1,4);
for i = 1:4
    [img, ~, alpha]  = imread(strcat(prm.stim.path, num2str(selected_img(i)), '.png')); 
    img              = imresize(img,prm.stim.size);
    alpha            = imresize(alpha,prm.stim.size);
    img(:, :, 4)     = alpha;
    Texture_img(i)   = Screen('MakeTexture', ptb.win, img);
end
prm.stim.grating = img;

subscreen_rects = cell(1,4);
for icenter = 1:4
    rects = zeros(4, 4);
    for quadrant = 1:4    
        bias_x = prm.stim.subcenters(icenter,1) + prm.stim.disControl(quadrant) * prm.stim.distance(quadrant,1);
        bias_y = prm.stim.subcenters(icenter,2) + prm.stim.disControl(quadrant) * prm.stim.distance(quadrant,2);
        rects(quadrant,:) = MakeOffsetRect(ptb, prm.stim, bias_x, bias_y,quadrant, 0);
    end
    subscreen_rects{icenter} = rects;
end


%% background
prm.bg.path          = "/Users/songyunbai/Documents/Donders/Experiment/RIFT test/dataset/scrambled_bg_gray.jpg";
bg                   = imresize(imread(prm.bg.path), [ptb.win_w/2, ptb.win_h/2]); 
prm.bg.grating       = bg;
Texture_bg           = Screen('MakeTexture', ptb.win, bg);

% initialize drawing rects (centered on physical screen)
bgs   = zeros(4, 4);
for quadrant = 1:4
    bgs(quadrant,:) = MakeOffsetRect(ptb, prm.bg, -ptb.win_w/2, -ptb.win_h/2, quadrant,1);
end

%% fixation
% also store fixation dot as texture (easier drawing/translation later)
prm.stim.param.fix_rad_outer = 0.2;
prm.stim.param.fix_rad_inner = 0.1;
stim = prm.stim;
siz_px = ptb.win_w/2;
% draw the fixation dot into an anti-aliased buffer first, then copy into
% regular window (necessary to get anti-aliasing)
buf = Screen('OpenOffscreenWindow', ptb.win, [0 0 0 0], [0 0 siz_px siz_px], [], [], 4);
rect = [siz_px/2-stim.deg2pix(stim.param.fix_rad_outer), siz_px/2-stim.deg2pix(stim.param.fix_rad_outer), siz_px/2+stim.deg2pix(stim.param.fix_rad_outer), siz_px/2+stim.deg2pix(stim.param.fix_rad_outer)];
Screen('FillOval', buf, 255, rect);
rect = [siz_px/2-stim.deg2pix(stim.param.fix_rad_inner), siz_px/2-stim.deg2pix(stim.param.fix_rad_inner), siz_px/2+stim.deg2pix(stim.param.fix_rad_inner), siz_px/2+stim.deg2pix(stim.param.fix_rad_inner)];
Screen('FillOval', buf, 0, rect);
stim.tex.fix = Screen('OpenOffscreenWindow', ptb.win, [0 0 0 0], [0 0 siz_px siz_px]);
Screen('CopyWindow', buf, stim.tex.fix);

%% diode tracking square
prm.stim.param.diode_track_size = 3; % vis deg
sz = round(stim.deg2pix(prm.stim.param.diode_track_size));
stim.tex.diode_track = Screen('OpenOffscreenWindow', ptb.win, [0 0 0 0], [0 0 siz_px siz_px]);
Screen('FillRect', stim.tex.diode_track, 255, [0 siz_px-sz sz siz_px]);

%% ========================================================================
% tagged presentation loop
vbl = []; % the flip vbl timestamp
fliptimes = nan(1, size(prm.stim.tag_sigs,2)/12);
flipind = 1;
% tag_sig = repmat([0 0.25 0.5 1], [1 6]);
for phys_frame_ind = 1:size(prm.stim.tag_sigs,2)
    % where are we going to draw it?
    % quadrant will increase by 1 every 3 physical frames and reset after
    % 12 physical frames
    quadrant = mod(floor((phys_frame_ind-1)/3), 4) + 1;
    
    % select the proper colour channel to draw into
    % colorchan will increase by 1 every physical frame and reset after 3
    % physical frames
    colorchan = mod(phys_frame_ind-1, 3) + 1;
    colmask = zeros(4,1);
    colmask(colorchan) = 1;
    Screen('BlendFunction', ptb.win, [], [], colmask);
    
    % draw the stimulus with the specified tagging type

    % luminance tagging between 0 and 100%: multiply intensity with tagging signal
    % (black stays black, white goes between 0-100%)
    
    Screen('DrawTexture', ptb.win, Texture_bg,[],bgs(quadrant,:));
    Screen('DrawTextures', ptb.win, stim.tex.fix, [], bgs');
    for isti = 1:4
        color_mod = ones(1, 3) * 255 * prm.stim.tag_sigs(isti, phys_frame_ind);
        sti_pos = subscreen_rects{quadrant};
        Screen('DrawTexture', ptb.win, Texture_img(isti),[],sti_pos(isti,:),[],[],[], color_mod);
        Screen('DrawTexture', ptb.win, stim.tex.diode_track, [],...
        sti_pos(isti,:), [], [], [], 255*prm.stim.tag_sigs(isti, phys_frame_ind));
    end

%     Screen('DrawTexture', ptb.win, Texture_img(quadrant),[],rects(quadrant,:),[],[],[], color_mod);
  
        
    % draw diode tracking square
%     Screen('DrawTexture', ptb.win, stim.tex.diode_track, [],...
%         diode_rects(quadrant,:), [], [], [], 255*tag_sig(phys_frame_ind));
    
    % flip if necessary (every 12 physical frames)
    if phys_frame_ind > 1 && mod(phys_frame_ind, 12) == 0
        % draw fixation dot in every quadrant
%         Screen('BlendFunction', ptb.win, [], [], mask_old);
%         Screen('DrawTextures', ptb.win, stim.tex.fix, [], rects');
        if ~isempty(vbl)
            vbl = Screen('Flip', ptb.win, vbl+ptb.ifi/2);
        else
            % the first time we're flipping we don't have a previous
            % timestamp yet
            vbl = Screen('Flip', ptb.win);
        end
        fliptimes(flipind) = vbl;
        flipind = flipind + 1;
    end

%     % for debugging 
%     doclear = phys_frame_ind > 1 && mod(phys_frame_ind, 12) == 0;
%     fprintf('flipping, tag=%.4g %%...', tag_sig(phys_frame_ind)*100);
%     Screen('Flip', ptb.win, 0, double(~doclear));
%     fprintf('flipped.\n');
%     
%     % store image
%     allimgs{end+1} = Screen('GetImage', ptb.win);
%     WaitSecs(0.1);
% %     KbWait();
end

for i = 1:4
    Screen('Close', Texture_img(i)) 
end
Screen('Close', Texture_bg) 
Screen('Close', stim.tex.fix) 
Screen('Close', stim.tex.diode_track) 
sca
