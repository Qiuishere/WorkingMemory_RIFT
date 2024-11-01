function monitor = open_screen(monitor,screen_size,background,gamma,translate,scale)
% wrapper to Psychtoolbox Screen
%                                                     D.Pascucci 09.11.2021
%--------------------------------------------------------------------------
default('screen_size',monitor.screen_pix)
default('background',monitor.grey)
default('gamma',[])
default('translate',[])
default('scale',[])
% open window
[monitor.window, ~]    = PsychImaging('OpenWindow', monitor.screen_num,...
    background,screen_size,...
    [],[],[],[],4); % multisampling
if ~isempty(gamma)
    Screen('LoadNormalizedGammaTable',monitor.window,repmat(gamma,[1 3]));
end
monitor.ifi            = Screen('GetFlipInterval', monitor.window);
Screen('BlendFunction', monitor.window, 'GL_SRC_ALPHA', ....
    'GL_ONE_MINUS_SRC_ALPHA');

if ~isempty(translate)
    Screen('glTranslate', monitor.window, translate(1)*monitor.width, translate(2)*monitor.height);
end
if ~isempty(scale)
    Screen('glScale', monitor.window, scale(1), scale(2));
end

% other defaults here
Screen('TextSize', monitor.window, 20);
HideCursor;
Priority(max(Priority));

% % update fields when the screen size changes (NOW IN monitor_info.m)
% if nnz(tmp~=monitor.screen_pix)>0
%     original           = monitor;
%     % overwrite
%     [monitor.width, monitor.height] = deal(screen_size(3),screen_size(4));
%     monitor.original   = original;
%     monitor.rate       = 1/monitor.hz;
%     monitor.frame_ms   = @(ms) round(ms/(1000/monitor.hz));
%     [width, height]    = Screen('DisplaySize', original.screen_num);
%     monitor.screen_cm  = .1*[width, height];
%     monitor.screen_pix = [0 0 monitor.width monitor.height];
%     [x, y]             = RectCenter(monitor.screen_pix);
%     monitor.center     = [x, y];
%     monitor.screen_num = monitor.original.screen_num;
%     monitor.window     = monitor.original.window;
%     % white, black and gray
%     monitor.black      = BlackIndex(monitor.screen_num);
%     monitor.white      = WhiteIndex(monitor.screen_num);
%     monitor.grey       = monitor.white / 2;
%     % pixels, cm and degrees
%     monitor.pix_in_cm  = monitor.screen_cm(1)/monitor.width;
%     monitor.deg_to_cm  = 2*monitor.original.view_dist*tan(pi*1/(2*180));
%     monitor.deg_to_pix = fix(monitor.deg_to_cm/monitor.pix_in_cm);
% end

% TO CHECK (substituting line 5)
% [monitor.window, tmp]  = PsychImaging('OpenWindow', monitor.screen_num, ...
%                                       background,screen_size,32,2,    ...
%                                       [],[],kPsychNeed32BPCFloat);
% [monitor.window , tmp] = Screen('OpenWindow', monitor.screen_num,...
%                                       background, screen_size);
