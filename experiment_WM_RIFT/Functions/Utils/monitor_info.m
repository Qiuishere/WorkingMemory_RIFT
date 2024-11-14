function monitor = monitor_info(screenId,Environment)

%%remember to change the parameters to the correct ones in the lab
default('screenId',max(Screen('Screens')))

% initialize
PsychDefaultSetup(2);

% get screen into and confirm it's in the correct setting
monitor            = Screen('resolution',  screenId);
monitor.screenId   = screenId;

switch Environment
    
    case 1 % Office        
        correctRefRate = 60;
        monitor.viewDist = 120;
        
    case 2 % In MEG
        correctRefRate = 120;
        correctRes = [1280 720]; % This is the Resolution we are supposed to use
        monitor.viewDist = 100; % viewing distance in cm
        
    case 3 % In Behaviour lab
        correctRefRate = 120;
        correctRes = [1280 720]; % This is the Resolution we are supposed to use
        monitor.viewDist = 95; % viewing distance in cm
   
end

if ismember(Environment,[2,3])
    if (monitor.hz - correctRefRate) >= 1
        error('Framerate is not %g!', correctRefRate);
    end
    if ~isequal([monitor.width,monitor.height], correctRes)
        error('Resolution is not %g x %g!', correctRes(1), correctRes(2))
    end
end

monitor.dt         = 1/monitor.hz;
monitor.frame_ms   = @(ms) round(ms/(1000/monitor.hz));
[width, height]    = Screen('DisplaySize', screenId);

monitor.screen_cm  = .1*[width, height];
monitor.screen_pix = [0 0 monitor.width monitor.height];
[x, y]             = RectCenter(monitor.screen_pix);
monitor.center     = [x, y];
monitor.screen_num = screenId;
% white, black and gray     
monitor.black      = BlackIndex(screenId);
monitor.white      = WhiteIndex(screenId);
monitor.grey       = monitor.white / 2;
% pixels, cm and degrees
monitor.pix_in_cm  = monitor.screen_cm(1)/(monitor.width/2);
monitor.deg_to_cm  = 2*monitor.viewDist*tan(pi*1/(2*180));
monitor.pix_to_deg = monitor.screen_cm(1)/(monitor.width/2)/ ...
                     monitor.viewDist*180/pi;
monitor.deg_to_pix = fix(1/monitor.pix_to_deg);
% additional
info.struct        = Screen('Version');
info.comp          = Screen('Computer');
info.arch          = computer;
monitor.info       = info;
% =========================================================================
